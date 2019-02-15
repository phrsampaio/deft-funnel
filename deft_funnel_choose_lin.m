function [ sampleSet, iterate, evaluations, xstatus, sstatus ] =            ...
    deft_funnel_choose_lin( f, c, sampleSet, iterate, setting,              ...
    evaluations, modelSize, Delta, sstatus, const )
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Builds a safely nondegenerate set of points Y of degree n+1 using the 
% points from X (the set of all points) which are contained in the current 
% Delta and which lie in the current subspace. If the current space is the 
% full space then no dummy points from X are considered for Y. To choose 
% appropriate points from X, a greedy algorithm is used. If less then n+1 
% points were found using this procedure, the set Y is augmented by new 
% optimally placed points.
%
% If only one single point is available in X, the vertices of a simplex are 
% taken to construct Y.
%
% If any new points were computed, they are also added to X and their function
% values are computed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
n           = iterate.xdim;               
xstatus     = zeros( sampleSet.nbPoints, 1 ); % initialize xstatus (points from X in Y)
lb          = setting.lx( iterate.indfree );
ub          = setting.lx( iterate.indfree );
I           = eye( n );
Ifull       = eye( iterate.fulldim );
CNTsin      = setting.CNTsin + 1;

sampleSet.Y = []; % Sample set Y will be built

if ( setting.verbose > 1 )
    verbose = 1;
else
    verbose = 0;
end

if ( verbose )
    disp(' *** enter deft_funnel_choose_lin ***************************')
end

% Define subset of X inside Delta and in the current subspace
k = 0;
for i = 1:sampleSet.nbPoints
   if ( (norm(iterate.x - sampleSet.X( iterate.indfree, i )) <= ...
          Delta * (1 + setting.eps_L)) && (sstatus(i) == const.in) )

      k = k + 1;
      XinDelta(:,k)   = sampleSet.X( iterate.indfree, i );
      ind_XinDelta(k) = i;
      
      if ( i == sampleSet.i_xbest )
         i_xbestinDelta = k; 
      end
   end
end 

% Find a safely nondegenerate set of points
if ( length( ind_XinDelta ) > 1 )
   strategy                 = 'greedy';
   [ind_YinDelta, Y, pY ]   = deft_funnel_greedy( i_xbestinDelta, XinDelta, setting.kappa_th );
   sampleSet.Y              = Y;
   sampleSet.ind_Y          = ind_XinDelta(ind_YinDelta);
   xstatus(sampleSet.ind_Y) = const.inY;
else
   strategy                 = 'simplex';
   sampleSet.Y(:,1)         = XinDelta;
   sampleSet.ind_Y          = ind_XinDelta;
   xstatus(sampleSet.ind_Y) = const.inY;
   pY                       = 1;
end

if ( verbose )
   disp(' *** constructed safely nondegenerate set of points ***')
   strategy
end

% Return if n+1 points in the set
if ( pY >= setting.rep_degree )
   if (verbose)
      disp(' *** found n+1 good points in X ***')
      disp(' *** exit choose_lin *********************')
   end
   return
end

% Compute new points if degree less than required degree
kk = 0;
if ( pY < setting.rep_degree )

   augment = setting.rep_degree - pY;

   if ( strcmp( strategy, 'greedy' ) )

      while ( pY < setting.rep_degree )

         % Pick a quasi-random interpolation point
         kk = kk + 1;
         ynew = Delta*sin(kk*(CNTsin+[1:n])');

         [ sampleSet, pY ] = deft_funnel_augment_Y( sampleSet, ynew, setting );

         % Optimally replace it
         if ( setting.hardcons )
            [ ynew, ~, msgTR ] = deft_funnel_find_new_yj_bc( sampleSet,     ...
                iterate, pY, setting, Delta );
         else
            [ ynew, ~, msgTR ] = deft_funnel_find_new_yj( sampleSet, pY,    ...
                setting, Delta );
         end

         if ( strcmp( msgTR(1:6), 'Error0') )
            Y(:,pY) = [];
            pY = pY - 1; 
         else
            sampleSet = deft_funnel_replace_in_Y( sampleSet, ynew, pY, setting );
         end

      end

   elseif ( strcmp( strategy, 'simplex' ) )

      % Build a simplex if only one point in the given set
      I = eye( n );
      sampleSet.Y( :, 1 ) = iterate.x;
      
      % Check Delta for not being too big
      if ( setting.hardcons == 1 )
         for j=1:n
      
            temp = ub(j)-lb(j);
       
            if (temp < Delta + Delta)         
               Delta = 0.5 * temp;            
            end
                     
         end         
      end
      
      % Take alternating direction of the simplex vertices due to
      % diversity reasons using the variable CNTsin
      if ( mod(CNTsin,2) == 1 )
         step1 = Delta;
         step2 = -Delta;
      else
         step1 = -Delta;
         step2 = Delta;
      end

      % Take care that all points are inside the bounds
      for j = 1:n
         lb_out = 0; ub_out = 0;
      
         if ( setting.hardcons == 1 && x(j) + step1 < lb(j) )
            step1 = -step1;
            lb_out = 1;
         elseif ( setting.hardcons == 1 && x(j) + step1 > ub(j) )
            step1 = -step1;
            ub_out = 1;
         end
         
         sampleSet.Y( :, j+1 ) = iterate.x + step1 * I( :, j );

         if ( setting.rep_degree >= modelSize.pdiag )
         
            if ( setting.hardcons == 1 && ( x(j) + step2 < lb(j) || lb_out == 1 ) )
               step2 = min( 2*Delta, ub(j) - x(j) );
            elseif ( hardcons == 1 && ( x(j) + step2 > ub(j) || ub_out == 1 ) )
               step2 = max( -2*Delta, lb(j)-x(j) );
            end
            sampleSet.Y( :, j+1+n ) = x + step2 * I( :, j );
         end
         
      end
      
      sampleSet = deft_funnel_build_QR_of_Y( sampleSet, setting );
      
   end % end of strategies

   sample = deft_funnel_create_sample( [], [], iterate );
   
   for i = 1:augment
      
      k = setting.rep_degree - augment;

      % Set index of new point
      sampleSet.ind_Y(k+i) = sampleSet.nbPoints+1;
      sample.x = sampleSet.Y( :, k+i );

      % Update X and evaluate function at ynew
      [ sampleSet, sample, evaluations, xstatus, sstatus ] =        ...
       deft_funnel_augmX_evalfc( f, c, sample, sampleSet, setting,  ...
       evaluations, xstatus, const.inY, sstatus );
   end

end

sampleSet.i_xbest  = sampleSet.ind_Y( 1 );
sampleSet.fY       = sampleSet.fX( sampleSet.ind_Y );
sampleSet.cY       = sampleSet.cX( :, sampleSet.ind_Y );

if ( verbose )
   disp(' *** augmented Y and X with new point(s) ***')
   disp(' Y = ')
   sampleSet.Y
   disp(' *** exit deft_funnel_choose_lin *********************')
end

end % end of deft_funnel_choose_lin
