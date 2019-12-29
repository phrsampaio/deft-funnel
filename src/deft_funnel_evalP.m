function value = deft_funnel_evalP(P, x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Computes the value of the polynomial P at x, where P is represented by
% the row vector containing its coefficients for the successive monomials.
% More specifically, these values are:
% P(1)         : constant coefficient,
% P(2:n+1)     : coefficients for the linear terms in x(1)... x(n),
% P(n+2:2n+2)  : coefficients for the squared terms in x(1)^2 ... x(n)^2
% P(2n+3,3n+2) : coefficients for the quadratic terms of the first subdiagonal: 
%                in x(1)*x(2) ... x(n-1)*x(n)
% (3n+3,4n+1)  : coefficients for the quadratic terms of the second subdiagonal: 
%                in x(1)*x(3) ... x(n-2)*x(n)
% etc.
%
% Input:
%   - P        : the polynomial model
%   - x        : the point at which the model must be evaluated
%
% Output:
%   - value    : the value of the model at x.
%
% Dependencies : deft_funnel_evalZ
% Called by    : deft_funnel_evalL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

value = P * deft_funnel_evalZ(x, size(P, 2));

