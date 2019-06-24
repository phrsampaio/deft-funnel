function [ sampleSet, replaced, maximprove, badcond ] = ...
           deft_funnel_repair_Y( sampleSet, iterate, setting, Delta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Repairs the interpolation set Y by first replacing the interpolation points
% at a distance larger than setting.factor_FPR*Delta from Y(:,1), and then replaces 
% interpolation points with the best possible point obtained by taking the 
% best optimal replacement, until this improvement no longer causes a relative 
% simplex volume increase by at least threshold.  For conventions on how 
% polynomials are represented, see the documentation of evalZ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_improve_loops = 20;        % the max number of times every close point is 
                               % considered for improvement
badcond           = 0;         % interpolation matrix conditioning
                               % 0 = well conditioned; 1 = ill conditioned

if ( setting.verbose > 1 )
    verbose = 1;
else
    verbose = 0;
end

if ( verbose > 0 )
    disp([ ' --- enter deft_funnel_repair_Y for Delta = ', num2str(Delta), ' ---' ] );
    disp(' ')
end

[ ~, p1 ] = size( sampleSet.Y );
replaced  = [];

% Compute the distance of the current interpolation points to the base point.
d = zeros( 1, p1 );
for j=2:p1
   d(j) = norm( sampleSet.Y(:,j)- sampleSet.Y(:,1) );
end
[ dsorted, jsorted ] = sort( d, 'descend' );

if ( verbose )
   disp(' Distance of current interpolation points to the iterate: ')
   d
   dsorted
   jsorted
end

% First replace the distant points.
if ( verbose )
   disp(' Trying to replace far points first ')
   disp(' ')
end

jmax = 0;
for j = 1:p1
   if ( dsorted(j) > setting.factor_FPR*(1 + setting.eps_L)*Delta )
      jmax = jsorted(j);
      
      if ( setting.hardcons == 1 )
         [ y, improvement ] = deft_funnel_find_new_yj_bc( sampleSet,        ...
             iterate, jmax, setting, Delta );
      else
         [ y, improvement ] = deft_funnel_find_new_yj( sampleSet, jmax,     ...
             setting, Delta );
      end
      
      if ( verbose > 0 )
        disp([' ==> lambda = ' num2str(improvement) ' at j=' num2str(jmax)])
      end
      
      if ( improvement > setting.Lambda_FP )
         
         sampleSet = deft_funnel_replace_in_Y( sampleSet, y, jmax, setting );
         replaced  = [ replaced jmax ];
         d( jmax ) = norm( y - sampleSet.Y(:,1) ); 
      end
      
   else
      sampleSet.Y_radius = dsorted(j);
      break
   end
   
end

if ( verbose > 0 )
    oldLambda = sampleSet.lambda;
    oldRadius = sampleSet.Y_radius;
	sampleSet = deft_funnel_poisedness_Y( sampleSet, iterate, setting );
	disp(' Results of attempt to replace far points:')
    disp([' ==> Number of points replaced: ', num2str(length(replaced))] )
    replaced
    disp([' ==> poisedness(Y) before replacement = ', num2str(oldLambda)])
    disp([' ==> poisedness(Y) after replacement = ', num2str(sampleSet.lambda)])
    disp([' ==> sampleSet.Y_radius before replacement  = ', num2str(oldRadius) ])
    disp([' ==> sampleSet.Y_radius after replacement  = ', num2str(sampleSet.Y_radius) ])
    disp(' ')
end

% Perform a loop over possible optimal improvements.
if ( verbose )
   disp(' Trying to replace any point based on the maximum possible improvement')
   disp(' ')
end

for  k = 1:max_improve_loops

   maximprove = 0;

   % Loop on all the possible replacements to find the best one.
   for j = 2:p1

      if ( setting.hardcons == 1 )
         [ y, improvement ] = deft_funnel_find_new_yj_bc( sampleSet,        ...
             iterate, j, setting, Delta );
      else
         [ y, improvement ] = deft_funnel_find_new_yj( sampleSet, j,        ...
             setting, Delta );
      end

      % Remember the current polynomial value, index and replacement point
      % that is the best so far
      if ( verbose > 0 )
         disp([ ' ==> j = ', int2str(j), ' improve = ', num2str(improvement)])
         y
      end
      
      if ( improvement > maximprove )
         maximprove = improvement;
         jmax       = j;
         ymax       = y;
      end
   end

   % If no significant improvement was found, return after updating the 
   % interpolation radius.
   if ( maximprove < setting.Lambda_CP || jmax == 0 )
      sampleSet.Y_radius = max( d );
      if ( verbose > 0 )
         disp(' No significant improvement was found')
         disp([' maximprove(small) = ', num2str(maximprove), ...
               ', jmax= ' int2str(jmax)])
      end
      return
   end

   % Perform the best replacement.
   if ( verbose > 0 )
      disp([' maximprove= ' num2str(maximprove) ', jmax= ' int2str(jmax)])
   end

   sampleSet = deft_funnel_replace_in_Y( sampleSet, ymax, jmax, setting );
   d( jmax ) = norm( ymax - sampleSet.Y(:,1) );

   if ( length( find( replaced == jmax ) ) == 0 )
      replaced = [ replaced jmax ];
   end

end

% Recompute sampleSet.Y_radius after the exchanges
sampleSet.Y_radius = max( d ); 

if ( verbose > 0 )
   sampleSet = deft_funnel_poisedness_Y( sampleSet, iterate, setting );
   disp(' Final results of repairing Y:')
   disp([' ==> Number of points replaced: ', num2str(length(replaced))] )
   replaced
   disp([' ==> poisedness(Y) = ', num2str(sampleSet.lambda)])
   disp([' ==> sampleSet.Y_radius  = ', num2str(sampleSet.Y_radius) ])
   disp(' --- exit deft_funnel_repair_Y ---')
end

end % end of deft_funnel_repair_Y