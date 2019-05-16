function [ sampleSet, iterate, evaluations, xstatus, sstatus, poised_model, ...
           msg ] = deft_funnel_build_initial_sample_set( f, c,              ...
           sampleSet, iterate, setting, evaluations, xstatus, sstatus,      ...
           const, modelSize, Delta )
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Builds the initial sample set using either of the followinbg 
% two stratagies:
% (1) random interpolation points;
% (2) simplex approach.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( setting.verbose >= 2 )
   disp( ' *********** deft_funnel_build_initial_sample_set ***********' )
   disp( ' *** Building initial sample set' )
   disp( [ ' *** Degree of the initial model = ', int2str( setting.cur_degree ) ] )
end

% Compute an initial poised interpolation set around the starting point
sampleSet.Y( :, 1 ) = iterate.x;
if ( strcmp( setting.initial_Y, 'random' ) )

   sampleSet.Y( :, 2:setting.cur_degree ) = -ones( iterate.xdim, setting.cur_degree - 1 ) +                 ...
       2*rand( iterate.xdim, setting.cur_degree - 1 );

   for  j = 2:setting.cur_degree
      sampleSet.Y( :, j ) = sampleSet.Y( :, 1 ) + ...
           sampleSet.Y( :, j ) * ( Delta / norm( sampleSet.Y( :, j ) ) );
      
      % Make sure that the interpolation points are within the bounds
      sampleSet.Y( :, j ) = max( min( sampleSet.Y( :, j ), setting.ux ), setting.lx );
   end
   
   % Build the initial factorization
   sampleSet = deft_funnel_build_QR_of_Y( sampleSet, setting );

   % Make the set poised
   sampleSet = deft_funnel_repair_Y( sampleSet, iterate, setting, Delta );

elseif strcmp( setting.initial_Y, 'simplex' )

   % Compute the initial interpolation set (simplex plus midpoints)
   I = eye( iterate.xdim );

   for j = 1:iterate.xdim

      % Linear
      step1 = -Delta;
      sampleSet.Y( :, j + 1 ) = iterate.x + step1 * I( :, j );

      % Diagonal
      if ( setting.cur_degree >= modelSize.pdiag )
         step2 =  Delta;
         sampleSet.Y( :, j + 1 + iterate.xdim ) = iterate.x + step2 * I( :, j );
      end
   end

   % Quadratic
   if ( setting.cur_degree == modelSize.pquad )
      k = 2 * iterate.xdim + 2;
      for j = 1:iterate.xdim
         for jj = j+1:iterate.xdim
            sampleSet.Y( :, k ) = 0.5 * ( sampleSet.Y( :, j + 1 ) + sampleSet.Y( :, jj + 1 ) );
            k = k + 1;
         end
      end
   end
   
   % Make sure that the interpolation points are within the bounds
   for  j = 2:setting.cur_degree
      sampleSet.Y( :, j ) = max( min( sampleSet.Y( :, j ), setting.ux ), setting.lx );
   end
   
   % Build the initial factorization
   sampleSet = deft_funnel_build_QR_of_Y( sampleSet, setting );
end 

% Compute the models' gradient error associated to the sample set by
% computing the poisedness of the sample set
sampleSet = deft_funnel_poisedness_Y( sampleSet, iterate, setting );

% Compute the associated function values, possibly reducing Delta to
% ensure that the objective function remains finite at all interpolation
% points.

% Set a sample structure
sample = deft_funnel_create_sample( [], [], iterate );

for i = 1:setting.cur_degree

   % Store the new points in X and evaluate f and c at them
   sample.x = sampleSet.Y( :, i );
   [ sampleSet, sample, evaluations, xstatus, sstatus ] =        ...
       deft_funnel_augmX_evalfc( f, c, sample, sampleSet, setting,        ...
       evaluations, xstatus, const.inY, sstatus );

   sampleSet.fY( i )    = sampleSet.fX( i );
   sampleSet.cY( :, i ) = sampleSet.cX( :, i );

   if ( evaluations.nfeval > setting.maxeval || evaluations.nceval > setting.maxeval )
       msg  = [' Error: Maximum number of ', int2str( setting.maxeval ), ' function evaluations reached.'];
       % Including fixed variables at return
       if ( iterate.nfix > 0 )
           I  = eye( iterate.xdim + iterate.nfix );
           iterate.x  = I( :, iterate.indfix ) * setting.lx( iterate.indfix ) + I( :, iterate.indfree ) * iterate.x;
       end
       disp(msg)
       poised_model = 0;
       return
   end

end

xstatus = xstatus';
poised_model = 1;

% Retrieve function values of the iterate
iterate.feval      = sampleSet.fY( 1 );
iterate.ceval      = sampleSet.cY( :, 1 );

sampleSet.ind_Y    = 1:setting.cur_degree;
sampleSet.i_xbest  = sampleSet.ind_Y( 1 );

msg = 'Initial sample set successfully built';

if ( setting.verbose >= 2 )
   disp( ' *** Initial sample set successfully built' )
   disp( ' ***** return from deft_funnel_build_initial_sample_set *****' )
end

end % end of deft_funnel_build_initial_sample_set