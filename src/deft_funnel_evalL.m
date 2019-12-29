function values = deft_funnel_evalL(sample_set, x, setting, choice_set)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Computes the values of the Lagrange polynomials at x.
%
% Input:
%   - sample_set : struct of the sample set
%   - x          : the point at which the polynomial must be evaluated
%   - setting    : struct of parameters
%   - choice_set : indices (of Y's columns) for which to calculate the 
%                  Lagrange polynomial values
%
% Output:
%   - values     : the values of the Lagrange polynomials at x.
%
% Dependencies   : deft_funnel_evalP
% Called by      : deft_funnel_include_in_Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n, p1] = size(sample_set.Y);
I       = eye(p1);
lc      = length(choice_set);
q       = ((n+1)*(n+2))/2;

if (setting.whichmodel == 3 && p1 < q)

   % For underdetermined regression model, use min l2-norm model
   setting.whichmodel = 2;

end

if (setting.whichmodel == 0)
    
   % Evaluate (sub-basis) Lagrange polynomial (p1 = q)
   warning off
   values(choice_set) = deft_funnel_evalP(I(choice_set,:) * (sample_set.RZ \ sample_set.QZ'), x);

elseif (setting.whichmodel == 1)

   if (p1 == n+1 || p1 == q)

   % Evaluate Minimum L2 norm Lagrange polynomial (p1 == n+1 or p1 == q)
      warning off
      values(choice_set) = deft_funnel_evalP(I(choice_set,:) * (sample_set.QZ / sample_set.RZ'), x);

   else

      % Evaluate Minimum Frobenius norm Lagrange polynomial (p1 <= q)
      M  = bcdfo_evalZ(sample_set.Y, q)';
      MQ = M(:, n+2:q);
      phi = deft_funnel_evalZ(x, q);

      warning off
      values(choice_set) =  [I(choice_set,:) zeros(lc,n+1)] * (sample_set.QZ * (sample_set.RZ' \ [MQ * phi(n+2:q); phi(1:n+1)]));

   end

elseif (setting.whichmodel == 2) 
   
   % Evaluate Minimum L2 norm Lagrange polynomial (p1 <= q)
   if (p1 < q)

      warning off
      values(choice_set) = deft_funnel_evalP(I(choice_set,:) * (pinv(sample_set.RZ) * sample_set.QZ'), x);

   else

      warning off
      values(choice_set) = deft_funnel_evalP(I(choice_set,:) * (sample_set.RZ \ sample_set.QZ'), x);

   end

elseif (setting.whichmodel == 3)
   
   % Evaluate Regression Lagrange polynomial (p1 >= q)
   warning off
   values(choice_set) = deft_funnel_evalP(I(choice_set,:) * (sample_set.QZ * pinv(sample_set.RZ')), x);
   
end

end % end of deft_funnel_evalL
