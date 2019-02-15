function H = deft_funnel_hessP( P, x )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Computes the Hessian of the polynomial P at x, where P is represented by
% the row vector containing its coefficients for the successive monomials.
% More specifically, these values are:
% P(1)        : constant coefficient,
% P(2:n+1)    : coefficients for the linear terms in x(1)... x(n),
% P(n+2:2n+2) : coefficients for the squared terms in x(1)^2 ... x(n)^2
% P(2n+3,3n+2): coefficients for the quadratic terms of the first subdiagonal: 
%               in x(1)*x(2) ... x(n-1)*x(n)
% (3n+3,4n+1): coefficients for the quadratic terms of the second subdiagonal: 
%               in x(1)*x(3) ... x(n-2)*x(n)
% etc.
% 
% Input:
%   - P       : a row vector containing the coefficient of the polynomial
%   - x       : the vector at which the Hessian of P is to be evaluated (not really
%               needed except for finding the space dimension, since we consider at most
%               quadratic polynomials).
%
% Output:
%   - H       : the Hessian matrix of P at x
%
% Programming: Ph. Toint,S. Gratton and A. Troeltzsch, April 2009.
% (This version 22 VI 2009)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n     = length( x );
p1    = length( P );
nquad = p1 - n - 1;

if( nquad > 0 )

    % Diagonal
    ndiag = min( nquad, n );
    H     = diag( [ P( n+2:n+1+ndiag ) zeros( 1, n-ndiag )] );
    nquad = nquad - ndiag;

    % Subdiagonals
    if ( nquad > 0 )
       k = 2*n+1;
       for i = 1:n-1
          nsd = min( n-i, nquad);
          if ( nsd > 0 )
              for j = 1:nsd
                H( i+j, j ) = P( k+j ); 
                H( j, i+j ) = P( k+j ); 
              end
              k = k + nsd;
              nquad = nquad - nsd;
          end
          if ( nquad == 0 )
             break;
          end
       end
    end
else
    H = zeros( n, n );
end

end % end of deft_funnel_hessP