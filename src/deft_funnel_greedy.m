function [ ind_Y, Y, p, msg ] = deft_funnel_greedy( i_xbest, X, kappa )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: An implementation of greedy algorithm to get a maximally linearly 
% independent interpolation set based on the algorithm for the optimal basis 
% problem by O. Burdakov.
%
% Input: 
%   - i_xbest  : index of starting point for the algorithm
%   - X        : set of available points
%   - kappa    : threshold for defining degeneracy of the set Y
%
%  Output:
%   - ind_Y    : indices of selected points in Y
%   - Y        : the safely nondegenerate set of points
%   - p        : cardinality of the set Y
%   - msg      : (error) message
%
% Programming: A. Troeltzsch, November 2009.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
verbose  = 0;
msg      = '';
n        = size(X,1);       % problem dimension 
m        = size(X,2);       % nbr of points in X
p        = 0;               % nbr of points in Y
ind_XwoY = 1:m;             % indices in set XwoY = X\Y
A        = [];
Hadamard = 1;
G        = 1;
k        = 1;               % iteration counter

if ( verbose > 1 )
    disp('*** enter greedy *********************')
end

% Compute the distances between the points (in Delta)
X2        = X' * X;
D         = diag(X2);
distances = sqrt( D * ones(1,m) + ones(m,1) * D' - 2 * X2 );

% Set very small distances to Infinity to not consider them
distances( find( distances <= 1.0e-6 ) ) = Inf;

% Include index of current iterate x in Y and exclude it from X
ind_Y             = i_xbest;
p                 = 1;
ind_XwoY(i_xbest) = [];

% Initial printout
if ( verbose )
    disp( '      It       cond     ind   p      m' )
    fprintf( '      %2d  %.7e  %2d  %2d  %5d\n', 0, Hadamard, i_xbest, p, m )
end

% Try to include points in Y
while (( length(ind_Y) < n+1 ) & ( k < m ))

   A_old = A;
   G_old = G;

   % Select smallest distance between the point(s) in Y and XwoY
   for i=1:p
      nr = ind_Y(i);
      [mdist(i),mind_in_XwoY(i)] = min( distances(nr,ind_XwoY) );
      mind_in_X(i) = ind_XwoY(mind_in_XwoY(i));
   end
   [mindist,minind_in_mdist] = min(mdist);
   minind_in_X = mind_in_X(minind_in_mdist);

   % Store the two points with the smallest distance
   xi    = X(:,ind_Y(minind_in_mdist));
   xmp1  = X(:,minind_in_X);
   y     = xmp1 - xi;

   % Update Hadamard condition number by using update-formula
   if ( p >= 2 )

      [Q,R]    = qr(A_old,0);
      xProjmp1 = y - (Q*Q')*y;

      G = G * norm(xProjmp1)^2 / norm(y)^2;

      if ( G ~= 0 )
         Hadamard = 1/G;
      else
         Hadamard = 1e+30;
      end

   end

   % Checking linear independence after including the point xmp1
   if ( Hadamard < kappa )

      % Update A
      A     = [ A y / norm(y) ];

      % Add respective index to Y
      ind_Y = [ind_Y minind_in_X];
      p     = p + 1;
   else

      % Roll back iteration if current point shouldn't be included
      if ( verbose > 1 )
         disp('*** roll back iteration ***')
      end

      G = G_old;
   end

   % Update indices in set X\Y (remove index of xmp1)
   ind_XwoY(mind_in_XwoY(minind_in_mdist)) = [];

   % Iteration printout
   if ( verbose )
      fprintf( '      %2d  %.7e  %2d  %2d  %5d\n', k, Hadamard, minind_in_X, p, m )
   end

   k = k + 1;
end

if ( verbose > 1 )
   disp('*** return from greedy ***************')
end

Y = X(:,ind_Y);

end % end of deft_funnel_greedy