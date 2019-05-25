function [ nit, sampleSet, iterate, setting, indicators, evaluations,       ...
           models,  xstatus, sstatus, vstatus, sspace_save, xspace_save,    ...
           Delta, Delta_f, Delta_c, vmax, msg ] =                           ...
            deft_funnel_ident_active_bnds( f, c, nit, sampleSet, iterate,   ...
            setting, indicators, evaluations, models, xstatus, sstatus,     ...
            vstatus, sspace_save, xspace_save, const, modelSize, Delta,     ...
            Delta_f, Delta_c, Deltamax, vmax, poised_model, level )
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Identifies active and nearly-active bounds, builds the models in the 
% subspace of the remaining free variables and call deft_funnel_main 
% recursively.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
msg       = '';                              % return message
indfree   = find( vstatus == const.free );   % currently free variables
indfix    = find( vstatus >= const.fixed );  % currently fixed variables
n         = length(indfree);                 % nb. of free variables
nfix      = iterate.fulldim - n;             % nb. of fixed variables
nbcons    = iterate.sdim;                    % number of constraints
lx_free   = setting.lx( iterate.indfree );   % lower bounds in the current subspace
ux_free   = setting.ux( iterate.indfree );   % upper bounds in the current subspace

if ( setting.verbose > 1 )
    verbose = 1;
else
    verbose = 0;
end

if ( verbose > 0 )
    disp( ' ' )
    disp( ' ********* deft_funnel_ident_active_bnds *********' )  
    disp( [ ' *** nfree = ', num2str( n )] );
    disp( ' *** Searching for new active bounds' )
end

% Initiate new subspace
xnew      = iterate.x;
xbd       = zeros( n, 1 );
i_free    = []; % Bounds not active
I_L       = []; % Active lower bounds
I_U       = []; % Active upper bounds
nbr_ss    = size( sspace_save, 2 );

if ( strcmp( level, 'toplevel' ) )
    % Define a new subspace from scratch
    sspace = zeros( iterate.fulldim, 1 );
else
    % Define the new subspace taking into account the current (last) one
    sspace = sspace_save( :, nbr_ss );
end

% Identify active and nearly-active bounds.
% If a bound is nearly-active, move current point to the bound
% and save it as xnew. Save value of the bound in xbd.

trial_vstatus = vstatus;
for i=1:n
    ii = iterate.indfree( i );
    if ( iterate.x( i ) - lx_free( i ) <= setting.eps_bnd )
        I_L           = [ I_L i ];
        xnew( i )     = lx_free( i );
        xbd( i )      = lx_free( i );
        trial_vstatus( ii ) = const.fixed;
        sspace( ii )  = -1;
    elseif ( ux_free( i ) - iterate.x( i ) <= setting.eps_bnd )
        I_U           = [ I_U i ];
        xnew( i )     = ux_free( i );
        xbd( i )      = ux_free( i );
        trial_vstatus( ii ) = const.fixed;
        sspace( ii )  = 1;
    else
        i_free        = [ i_free i ];
    end
end
I_LandU = [ I_L I_U ];

% Return if no active bounds
if ( isempty( I_LandU ) ) 
   msg = ' No active bounds';
   if ( verbose > 0 )
      disp( ' *** No active bounds were found' )
      disp( ' *** return from deft_funnel_ident_active_bnds ***' )
      disp( ' ' )
   end
   return
else
   if ( verbose > 0 )
      I_L
      I_U
   end
end 

% Return if new subspace was already tried before from this position

% xspace_save is a (3 x nbr_ss) matrix such that:
% - the first 2 lines contain the index of the starting and existing points in the new subspace
% - the 3rd line indicates if the starting point was already explored in
%   other subspaces

xspace_save( 3, : ) = 0; % no subspace has already been explored a priori
for i = 1:size( sspace_save, 2 )
   if ( norm(sspace - sspace_save( :, i ) ) == 0 &&                         ...
          (sampleSet.i_xbest == xspace_save( 1, i ) ||                      ...
           sampleSet.i_xbest == xspace_save( 2, i ) ) )
    
      msg = ' *** Subspace already explored';
      xspace_save( 3, i ) = 1;
      if ( verbose > 0 )
         disp( msg )
         disp( ' *** return from deft_funnel_ident_active_bnds ***' )
         disp( ' ' )
      end
      return
   end
end   

% Save actual subspace and index of current iterate
sspace_save                = [ sspace_save sspace ];
xspace_save( 1, nbr_ss+1 ) = sampleSet.i_xbest; % index of starting point in X
xspace_save( 2, nbr_ss+1 ) = 0;                 % index of exiting point in X (temporary)
  
% Check if x was moved to the bound 
% and evaluate f at xnew if necessary
if ( norm( iterate.x - xnew ) > setting.eps_bnd )

   % Check if new point coincides with another point
   coincides = 0;
   for i = 1:sampleSet.nbPoints
      if ( norm(xnew - sampleSet.X( iterate.indfree, i ) ) < setting.eps_bnd )
         coincides = i;
         break
      end
   end
   if ( coincides )
      msg = [ ' *** Current iterate moved to the bounds ',                  ...
                 'but this point already exists' ]; 
      if ( verbose > 0 )
         disp( msg )
         xnew
         disp( ' *** return from deft_funnel_ident_active_bnds ***' )
         disp( ' ' )
      end
      return
   end

   % Update X and evaluate f and c
   xnewSample = deft_funnel_create_sample( xnew, [], iterate );
   
   [ sampleSet, xnewSample, evaluations, xstatus, sstatus ] =               ...
       deft_funnel_augmX_evalfc( f, c, xnewSample, sampleSet,               ...
       setting, evaluations, xstatus, 0, sstatus);

   % Check if new function value is smaller
   % and adjust indices if necessary
   if ( xnewSample.feval < sampleSet.fX( sampleSet.i_xbest ) )      
      if ( verbose > 0 )
         disp( [' *** i_xbest: ' num2str( sampleSet.i_xbest ) ] )
         disp( ' *** after moving best point to the bounds: ' )
         disp( [' *** i_xbest: ' num2str( sampleSet.nbPoints ) ])
      end
      xstatus( sampleSet.nbPoints ) = const.inY;
      sampleSet.i_xbest             = sampleSet.nbPoints;
      sampleSet.ind_Y( 1 )          = sampleSet.nbPoints;
      sampleSet.fY( 1 )             = sampleSet.fX( end );
      sampleSet.cY( :, 1 )          = sampleSet.cX( :, end );
      iterate                       = xnewSample;
      
   else
      msg = ' *** Current iterate moved to the bound but new fvalue not smaller';
      if ( verbose > 0 )
          disp( msg )
          disp( ' *** return from deft_funnel_ident_active_bnds ***' )
          disp( ' ' )
      end
      return
   end
   
end
   
% Return if all bounds active
if ( length( I_LandU ) == n ) 

   msg = ' *** All bounds are active';
   
   if ( verbose > 0 )
      disp( msg )
      disp( ' *** return from deft_funnel_ident_active_bnds ***' )
      disp( ' ' )
   end
   return
   
% There are some non-active bounds
else
    
   iterate.indfree = indfree;
   iterate.indfix  = indfix;
   iterate.xdim    = n;        
   iterate.nfix    = nfix;    
   vstatus         = trial_vstatus;
    
   % Identify points in X which lie in the current subspace
   
   % sstatus(i) is set to 1 for X(:,i) when evaluating the function at
   % X(:,i) in 'deft_funnel_augmX_evalfc'
   
   for j = 1:sampleSet.nbPoints
   
      % exclude far points from the procedure
      if ( norm( iterate.x - sampleSet.X( iterate.indfree, j ) ) > Delta )
         sstatus( j ) = const.out;
      else
          
         for i = 1:n
            ii = iterate.indfree( i );

            if ( find( I_L == i ) >= 1 )

               % break if variable is NOT close to active lower bound
               if ( abs( sampleSet.X( ii, j ) - setting.lx( ii ) ) > setting.eps_bnd )
                  sstatus( j ) = const.out;
               end

            elseif ( find( I_U == i ) >= 1 )

               % break if variable is NOT close to active upper bound
               if ( abs( setting.ux( ii ) - sampleSet.X( ii, j ) ) > setting.eps_bnd )
                  sstatus( j ) = const.out;
               end
            end
         end     
      end
   end

   % Collect indices of points in the subspace
   ind_inSubspc = find( sstatus >= const.in ); 
    
   if ( verbose > 0 )
     disp( [ ' *** Nb. of points at an active bound: ',                     ...
             num2str( length( ind_inSubspc ) ) ] )
     disp( ' *** Their indices: ')
     ind_inSubspc'
         
   end
    
   % Recod info from the current subspace before entering the new one
   indfree_old = iterate.indfree;
   indfix_old  = iterate.indfix;
   xfix_old    = iterate.xfix;
   
   % Merge all fixed variables and define new subspace
   I               = eye( iterate.fulldim );
   iterate.xfix    = I( :, iterate.indfix ) * iterate.xfix( iterate.indfix ) + I( :, iterate.indfree ) * xbd;
   iterate.indfree = find( vstatus == const.free );
   iterate.indfix  = find( vstatus >= const.fixed );
   iterate.x       = sampleSet.X( iterate.indfree, sampleSet.i_xbest );
   iterate.xdim    = length( iterate.indfree );
   iterate.nfix    = length( iterate.indfix );
   new_n           = iterate.xdim;

   % Adjust the constants related to the dimension of the problem
   rep_degree_before = setting.rep_degree;
    
   if ( setting.rep_degree == modelSize.plin )
      setting.rep_degree = new_n + 1; 
   elseif ( setting.rep_degree == modelSize.pdiag )
      setting.rep_degree = 2 * new_n + 1; 
   elseif ( setting.rep_degree == modelSize.pquad )
      setting.rep_degree = ( ( new_n + 1 ) * ( new_n + 2 ) ) / 2; 
   end

   modelSize.plin  = new_n + 1;
   modelSize.pdiag = 2 * new_n + 1;
   modelSize.pquad = ( ( new_n + 1 ) * ( new_n + 2 ) ) / 2;
 
   %  Construct interpolation set in the subspace (n+1 points)
   [ sampleSet, iterate, evaluations, xstatus, sstatus ] =                  ...
       deft_funnel_choose_lin( f, c, sampleSet, iterate, setting,           ...
       evaluations, modelSize, Delta, sstatus, const );
   
   sampleSet.fY = sampleSet.fX( sampleSet.ind_Y );
   sampleSet.cY = sampleSet.cX( :, sampleSet.ind_Y );
   setting.cur_degree = size( sampleSet.Y, 2 );
   
   % Build the initial factorization
   sampleSet = deft_funnel_build_QR_of_Y( sampleSet, setting );
   
   % Compute maximum gradient error estimation
   sampleSet = deft_funnel_poisedness_Y( sampleSet, iterate, setting );

   %  Compute necessary information for the new subspace
   [ models, iterate, indicators ] = deft_funnel_update_iteration( sampleSet, iterate, setting );
   
   if ( verbose > 0 )
      disp( ' ' )
      disp( [ ' *** New set of points in subspace (sample set Y) nfree = ', ...
              num2str( new_n )] ) 
      sampleSet.Y
      disp( ' *** Computed new model ' )
   end
   
   it_type = 'entering subspace';
   
   %  Iteration printout after computing the subspace model
   deft_funnel_printout( nit, evaluations, iterate, setting,                ...
    indicators, vmax, 0, 0, Delta, 0, it_type, msg, sampleSet.errg );

   if ( verbose > 0 )
      disp( ' ' )
      disp( [ ' *** Solve subspace problem (nfree = ', ...
              num2str(new_n) ') ********'])
      disp( ' ##################################################')
   end
   
   % Return if problem already converged (after evaluating f at possible 
   % new points)
   if ( indicators.norm_c_s <= setting.epsilon && indicators.pi_f <= setting.epsilon ) 
      if ( verbose > 0 )
         disp( ' *** Problem already converged' )
         disp( ' ' )
      end
   end 
         
   % Solve problem in the subspace
   [ exit_algo, nit, sampleSet, iterate, indicators, evaluations, models,   ...
    xstatus, sstatus, vstatus, sspace_save, xspace_save, const,             ...
    modelSize, Delta, Delta_f, Delta_c, vmax, msg ] =                       ...
    deft_funnel_main( f, c, nit-1, nit-1, sampleSet, iterate, setting,      ...
    indicators, evaluations, models, xstatus, sstatus, vstatus,             ...
    sspace_save, xspace_save, const, modelSize, Delta, Delta_f, Delta_c,    ...
    Deltamax, vmax, poised_model, 'subspace');

   if ( verbose > 0 )
      disp( ' ' )
      disp( [ ' *** Exit subspace problem (nfree = ',                       ...
              num2str(new_n) ') ********'])
      disp( ' ##################################################' )
      disp( ' ' )
   end
    
   % Store index of best point when exiting all subspaces and going to 
   % the toplevel (to avoid reentering the same subspace later on)
   nbr_ss_last = size( sspace_save, 2 );

   for j=nbr_ss+1:nbr_ss_last
	   xspace_save( 2, j ) = sampleSet.i_xbest;
   end
end

% Compute higher dimensional information
iterate.indfix  = indfix_old;
iterate.indfree = indfree_old;
iterate.xfix    = xfix_old;
iterate.x       = sampleSet.X( iterate.indfree , sampleSet.i_xbest );
iterate.xdim    = length( iterate.indfree );
iterate.feval   = sampleSet.fX( sampleSet.i_xbest );
iterate.ceval   = sampleSet.cX( :, sampleSet.i_xbest );

vstatus( iterate.indfree ) = const.free; % variable status

if ( strcmp( level, 'toplevel' ) )
    sstatus = ones( sampleSet.nbPoints, 1 );
end
n = length( iterate.indfree );

% Reset some constants
setting.rep_degree = rep_degree_before;
modelSize.plin     = n + 1;
modelSize.pdiag    = 2 * n + 1;
modelSize.pquad    = ( ( n + 1 ) * ( n + 2 ) ) / 2;

Delta = setting.epsilon;

% If maxeval is attained, return quickly
if ( evaluations.nfeval >= setting.maxeval )
    if ( strcmp( level, 'toplevel' ) )
       it_type = 'back to fullspace';
    else
       it_type = 'back to upper level';
    end
    deft_funnel_printout( nit, evaluations, iterate, setting, indicators,   ...
    vmax, 0, 0, Delta, 0, it_type, '', sampleSet.errg );
    return
end

[ sampleSet, iterate, evaluations, xstatus, sstatus ] =                     ...
    deft_funnel_choose_lin( f, c, sampleSet, iterate, setting,              ...
    evaluations, modelSize, Delta, sstatus, const );

if ( verbose > 0 )
  if ( strcmp( level, 'toplevel' ) )
     fprintf( '\n *** Build Y in the full space ********' )
  else
     fprintf( '\n *** Build Y in the higher space ********' )
  end
  sampleSet.Y
end

sampleSet.fY  = sampleSet.fX( sampleSet.ind_Y );
sampleSet.cY  = sampleSet.cX( :, sampleSet.ind_Y );
setting.cur_degree = size( sampleSet.Y, 2 );

% Build the initial factorization
sampleSet = deft_funnel_build_QR_of_Y( sampleSet, setting );

% Compute maximum gradient error estimation
sampleSet = deft_funnel_poisedness_Y( sampleSet, iterate, setting );

if ( verbose > 0 )
  fprintf(' *** Compute higher dimensional models ********\n\n')
end

% Compute necessary information for the new subspace
[ models, iterate, indicators ] = deft_funnel_update_iteration( sampleSet, iterate, setting );
   
if ( strcmp( level, 'toplevel' ) )
   it_type = 'back to fullspace';
else
   it_type = 'back to upper level';
end

% Iteration printout after computing the higher subspace model
deft_funnel_printout( nit, evaluations, iterate, setting, indicators,       ...
    vmax, 0, 0, Delta, 0, it_type, '', sampleSet.errg );

% Adjust Delta if not converged in the full-space
if ( ( indicators.norm_c_s > setting.epsilon || indicators.pi_f > setting.epsilon ) &&                         ...
      strcmp( level,'toplevel' ) )
  [ Delta_c, Delta_f, Delta ] = deal( max( 0.1 * indicators.pi_f, setting.epsilon ) );
end

msg = ' *** Identification of active bounds done successfully';

if ( verbose > 0 )
   disp(msg)
   disp( ' *** return of deft_funnel_ident_active_bnds **' )
   disp( ' ' )
   disp( ' ' )
end

end % end of deft_funnel_ident_active_bnds