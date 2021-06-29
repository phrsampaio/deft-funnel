function [xp, res, free, atlb, atub] =  deft_funnel_blls_spwmin(x, d, A, g, ...
    res, lb, ub, free, atlb, atub, verbose)
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Part of the BLLS solver. See main file 'deft_funnel_blls_exp.m'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Some algorithmic parameters
epsdzer = 1e-12;
epsfeas = 1e-12;
epsconv = 1e-10;

% Problem size
[m,n] = size(A);

% Initialization
inds      = [1:n];
latlb     = length(atlb);
latub     = length(atub);
lfree     = length(free);
xp        = x;
qp        = 0;
ak        = 0;

% Loop on the successive segments of the piecewise quadratic
for k = 0:n

% Compute the slope along the projected direction
   lk = g(free)'*d(free);

   % Print information on the new point, if requested
   if (verbose > 1)
      fprintf('       %4d             %+.4e  %+.4e %4d %4d %4d\n', ...
              k, qp, lk, lfree, latlb, latub)
   end

   % Check for termination is the slope is non-negative
   if (lk >= -epsconv)
      break
   end

   % Compute the unconstrained line minimizer along the projected
   % direction
   if ( k == 0 )
      Ad = A(1:m,free)*d(free);
   else
      if (length(natb) > 0)
         Ad = Ad - A(1:m,natb)*d(natb);
      end
   end
   Adsq = Ad'*Ad;
   ak   = - lk / Adsq;

   % Restrict the step the maximal one maintaining feasibility
   for kf = 1:length(free)
      kk = free(kf);
      if (d(kk) >= epsdzer)
         ak = min(ak, (ub(kk) - xp(kk)) / d(kk));
      elseif (d(kk) <= - epsdzer)
         ak = min( ak, (lb(kk) - xp(kk)) / d(kk));
      end
   end

   % Compute next iterate and update the gradient and the quadratic's value
   xp(free) = xp(free) + ak * d(free); 
   aAd      = ak * Ad;       
   g(free)  = g(free) + A(1:m,free)' * aAd;
   res      = res + aAd;
   qp       = qp + ak * lk + 0.5 * aAd' * aAd;
 
   % Determine which new bounds are active and add them to 
   % the current active set
   natlb     = free(find(abs(xp(free)-lb(free))<= epsfeas));
   natub     = free(find(abs(xp(free)-ub(free))<= epsfeas));
   natb      = [natlb natub];
   atlb      = [atlb natlb];
   atub      = [atub natub];
   atb       = [atlb atub];
   free      = inds;
   free(atb) = [];
   latlb     = length(atlb);
   latub     = length(atub);
   lfree     = length(free);

end % of the loop on the successive segments

end % end of deft_funnel_blls_spwmin