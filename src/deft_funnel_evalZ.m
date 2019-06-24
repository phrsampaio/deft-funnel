function Z = deft_funnel_evalZ( X, q )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Computes the matrix Z(X), where X is a matrix whose columns contains 
% points of the underlying space.  The vector Z(x) is, for a vector x, given 
% by the sequence of the values of the (at most quadratic) monomials taken at x.  
% More specifically, these values are:
% Z(x)(1)        : 1,
% Z(x)(2:n+1)    : the linear terms x(1)... x(n),
% Z(x)(n+2:2n+2) : the diagonal terms of the quadratic: x(1)^2 ... x(n)^2
% Z(x)(2n+3:3n+2): the first subdiagonal of the quadratic: 
%                  x(1)*x(2) ... x(n-1)*x(n)
% Z(x)(3n+3:4n+1): the second subdiagonal of the quadratic: 
%                  x(1)*x(3) ... x(n-2)*x(n)
% etc.
%
% Input:
%   - X          : the matrix whose columns contains the points at which the 
%                  monomials should be evaluated.
%   - q          : the number of monomials considered (q <= (n+1)*(n+2)/2)
%
% Output:
%   - Z          : the matrix Z(X), of size q x m.
%
% Programming: Ph. Toint, January 2009.
% (This version 12 IV 2010)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ n, m ] = size( X );             % [ dimension of the space, number of points in X ]
nlin     = min( n+1, q );         % number of constant and linear terms
nquad    = max( 0, q-nlin );      % number of quadratic terms
nlin     = nlin - 1;              % number of linear terms
Z        = zeros( q, m );

if ( q == 1 )
   Z = ones( 1, m );                       % constant terms
elseif ( q <= n+1 )
   Z = [ ones( 1, m ); X(1:nlin,1:m) ];    % constant and linear
else
   ndiag = min( n, nquad );
   Z     = [ ones( 1, m ); X(1:n,1:m); 0.5*X(1:ndiag,1:m).^2 ]; % same + diagonal
   nquad = nquad - ndiag;
   if ( nquad > 0 )
      for k = 1:n-1                        % the (i+1)-th subdiagonal
          nsd = min( n-k, nquad );
          if ( nsd > 0 )
             Z = [ Z; X(k+1:k+nsd,1:m).*X(1:nsd,1:m) ];
             nquad = nquad - nsd;
          end
          if ( nquad == 0 )
             break;
          end
      end
   end
end

end % end of deft_funnel_evalZ