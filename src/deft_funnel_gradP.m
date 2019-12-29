function g = deft_funnel_gradP(P, x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Computes the gradient of the polynomial P at x, where P is represented by
% the row vector containing its coefficients for the successive monomials.
% More specifically, these values are:
% P(1)        : constant coefficient,
% P(2:n+1)    : coefficients for the linear terms in x(1)... x(n),
% P(n+2:2n+2) : coefficients for the squared terms in x(1)^2 ... x(n)^2
% P(2n+3,3n+2): coefficients for the quadratic terms of the first subdiagonal: 
%               in x(1)*x(2) ... x(n-1)*x(n)
% (3n+3,4n+1) : coefficients for the quadratic terms of the second subdiagonal: 
%               in x(1)*x(3) ... x(n-2)*x(n)
% etc.
%
% Input:
%   - P       : a row vector contains the coefficients of the polynomial
%   - x       : the point at which the gradient must be evaluated
%
% Output:
%   - g       : the gradient of P at x
%
% Programming: Ph. Toint, S. Gratton and A. Troeltzsch, April 2009.
% (This version 30 IV 2009)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n  = length(x);
p1 = length(P);

% First-order terms in the polynomial
ng        = min(n, p1-1);
g         = zeros(n, 1);
g(1:ng) = P(2:ng+1);

% Second-order terms
nquad = p1 - n - 1;
if(nquad > 0)

    % Diagonal
    ndiag      = min(nquad, n);
    g(1:ndiag) = g(1:ndiag) + P(n+2:n+1+ndiag)'.*x(1:ndiag);
    nquad      = nquad - ndiag;

    % Subdiagonals
    if (nquad > 0)
       k = 2 * n + 1;
       for i = 1:n-1
          nsd = min(n-i, nquad);
          if (nsd > 0)
             g(i+1:i+nsd) = g(i+1:i+nsd) + P(k+1:k+nsd)' .* x(1:nsd);
             g(1:nsd)     = g(1:nsd)     + P(k+1:k+nsd)' .* x(i+1:i+nsd);
             k     = k + nsd;
             nquad = nquad - nsd;
          end
          if (nquad == 0)
             break;
          end
       end
    end
end

end % end of deft_funnel_gradP