function sample_set = deft_funnel_poisedness_Y(sample_set, iterate, setting)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Computes the poisedness of the interpolation set Y in a ball of radius
% Delta centered at Y(:,1). Poisedness is defined here as the maximum aboslute 
% value of the Lagrange polynomials taken over the ball and for all polynomials.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda    = 0;
[~, p1] = size(sample_set.Y);

if (setting.verbose >= 2)
        disp(' ' )
        disp(' *********** deft_funnel_poisedness_Y ***********')
        disp(' *** Computing the poisedness of the sample set')
end

% Compute the radius of the poisedness ball.
Y_radius = 0;
for j=2:p1
   Y_radius = max(Y_radius, norm(sample_set.Y(:,j) - sample_set.Y(:,1)));
end

% Loop on all the possible replacements to find the best one.
for j = 2:p1
    
   if (setting.hardcons == 1)
   
      [~, improvement] = deft_funnel_find_new_yj_bc(sample_set, iterate, j, ...
          setting, Y_radius);
               
   else

      [~, improvement] = deft_funnel_find_new_yj(sample_set, j, setting, Y_radius);
      
   end

   % Remember the current polynomial value, index and replacement point if
   % this is the best so far.
   lambda = max(improvement, lambda);

end

sample_set.lambda = lambda;
sample_set.Y_radius = Y_radius;
sample_set.errg = sample_set.lambda * sample_set.Y_radius;

if (setting.verbose >= 2)
        disp(' *** Poisedness of the sample set computed successfully')
        disp([' *** Gradient error estimation (sample_set.errg) =  ', int2str(sample_set.errg)])
        disp(' ***** return from deft_funnel_poisedness_Y *****')
        disp(' ')
end

end % end of deft_funnel_poisedness_Y