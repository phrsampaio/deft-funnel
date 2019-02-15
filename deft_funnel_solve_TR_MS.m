function [ s, lambda, norms, value, gplus, nfact, neigd, msg, hardcase ] =  ...
            deft_funnel_solve_TR_MS( g, H, Delta, eps_D )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: A simple implementation of exact trust-region minimization based on the
% Moré-Sorensen algorithm.
%
% Input: 
%   - g        : the model's gradient
%   - H        : the model's Hessian
%   - Delta    : the trust-region's radius
%   - eps_D    : the accuracy required on the equation ||s|| = Delta for a
%                boundary solution
% Output:
%   - s        : the trust-region step
%   - lambda   : the Lagrange multiplier corresponding to the trust-region constraint
%   - norms    : the norm of s
%   - value    : the value of the model at the optimal solution
%   - gplus    : the value of the model's gradient at the optimal solution
%   - nfact    : the number of Cholesky factorization required
%   - neigd    : the number of eifenvalue decompositions required
%   - msg      : an information message
%   - hardcase : 0 or 1
%
%  Dependencies: -
%
%  Programming: Ph. Toint and S. Gratton, April 2009. 
%  (This version 14 I 2010)
%
%  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose  = 0;
theta    = 1.0e-13;        % accuracy of the interval on lambda
epsilon  = 1.0e-15;        % theshold under which the gradient is considered as zero
nitmax   = 300;            % the maximum number of MS iterations
n        = length( g );    % space dimension
s        = zeros(n,1);     % initialize zero step
norms    = 0;              % initial stepnorm
lambda   = 0;              % initial lambda
value    = 0;              % initial model value
gplus    = g;              % initialize gplus
nfact    = 0;              % factorization counter
neigd    = 0;              % eigen decomposition counter
hardcase = 0;              % hard case

if ( verbose )
   disp( ' deft_funnel_solve_TR_MS : ============ enter' )
end

if ( length( find( isnan(H) ) ) ~= 0 )
   disp( ' deft_funnel_solve_TR_MS : H contains NaNs!')
   msg = 'error1';
   if ( verbose )
      disp( ' deft_funnel_solve_TR_MS : ============ error exit' )
   end
   return;
end

if ( length( find( ~isreal(H) ) ) ~= 0 )
   disp( ' deft_funnel_solve_TR_MS : H contains imaginary parts!')
   msg = 'error2';
   if ( verbose )
      disp( ' deft_funnel_solve_TR_MS : ============ error exit' )
   end
   return;
end

if ( length( find( isinf(H) ) ) ~= 0 )
   disp( ' deft_funnel_solve_TR_MS : H contains infinite elements!')
   msg = 'error3';
   if ( verbose )
      disp( ' deft_funnel_solve_TR_MS : ============ error exit' )
   end
   return;
end

% Compute initial bounds on lambda.
gnorm    = norm( g );
goverD   = gnorm / Delta;
Hnorminf = norm( H, inf );
if ( Hnorminf > 0 )        % because Octave generates a NaN for null matrices.
   HnormF   = norm( H, 'fro' );
else
   HnormF   = 0;
end   
lower    = max( 0, goverD - min( Hnorminf, HnormF ) );
upper    = max( 0, goverD + min( Hnorminf, HnormF ) );

% Compute the interval of acceptable step norms.
Dlower   = ( 1 - eps_D ) * Delta;
Dupper   = ( 1 + eps_D ) * Delta;

% Zero gradient
if ( norm( g ) < epsilon )
   if ( verbose )
      disp( ' deft_funnel_solve_TR_MS : ============ zero gradient:' )
   end
   [ V, D ]   = eig( H );
   neigd      = neigd + 1;
   [mu, imu ] = min( diag(D) );
   if ( mu < 0 )
      s = Delta * V(:,imu);
   else
       if ( norm(g) == 0 )
           s = zeros(size(g));
       else
           s = - Delta * ( g /norm (g) );
       end
   end
   sfound = 1;
   norms  = norm( s );
   lambda = -mu;

% Nonzero gradient
else

   if ( verbose )
      disp( ' deft_funnel_solve_TR_MS : ============ nonzero gradient:' )
   end

   % Compute initial lambda.
   if ( lower == 0 )
      lambda = 0;
   else
      lambda = max( sqrt( lower * upper ), lower + theta * ( upper - lower ) );
   end

   % Loop on successive trial values for lambda.
   for i = 1:nitmax
      new_lambda = -1;
      sfound     = 0;
      if ( verbose )
         disp( [ ' deft_funnel_solve_TR_MS (',int2str(i),'): lower = ',     ...
                 num2str( lower ), ' lambda = ', num2str( lambda ) ,        ...
                 ' upper = ',  num2str( upper ) ] )
      end

      % Factorize H + lambda I.
      [ R, p ] = chol( H + lambda * eye( n ) );
      if ( length( find( isnan( R ) ) ) ~= 0 )
         H
         lambda
         norm(g)
         R
         p
         disp( ' NaNs in Cholesky factorization' )
         msg = 'error4';
         if ( verbose )
            disp( ' deft_funnel_solve_TR_MS : ============ error exit' )
         end
         return;
      end
      nfact    = nfact + 1;

      % Successful factorization
      if ( p == 0 )
         warning off
         s      = - R \ ( R' \ g );
         sfound = 1;
         norms  = norm( s );
         if ( verbose )
            disp( [ ' deft_funnel_solve_TR_MS (',int2str(i),'): ||s|| = ', num2str( norms ), ...
                       ' Delta  = ', num2str( Delta ) ] )
         end

         % Test for successful termination.
         if ( ( lambda <= epsilon && norms <= Dupper ) || ...
              ( norms >= Dlower && norms <= Dupper ) )

            % Compute the optimal model value and its gradient.
            w     =  H * s;
            value = g' * s  + 0.5 * s' * w;
            gplus = g + w;
            norms = norm( s );
   
            % Define information message.
            if ( norms < ( 1 - eps_D ) * Delta )
               msg = 'interior solution';
            else
               msg = 'boundary solution';
            end
            if ( verbose )
               disp( ' deft_funnel_solve_TR_MS : ============ successful exit' )
            end
            return;
         end

         % Newton's iteration on the secular equation
         warning off
         w = R' \ s;
         normw2     = w'*w;
         new_lambda = lambda + ( ( norms - Delta ) / Delta ) * ( norms^2 / normw2 );

         % Check feasibility wrt lambda interval.
         if ( norms > Dupper )
            lower = lambda;
         else
            upper = lambda;
         end
         theta_range = theta * ( upper - lower );
         if ( new_lambda > lower + theta_range & new_lambda < upper - theta_range )
            lambda = new_lambda;
         else
            lambda = max( sqrt( lower * upper ), lower + theta_range );
         end

      % Unsuccessful factorization: take new lambda as the middle 
      % of the allowed interval
      else
         lower  = lambda;
         t      = 0.5;
         lambda = ( 1 - t ) * lower + t * upper;
      end

      % Terminate the loop if the lambda interval has shrunk to meaningless.
      if ( upper - lower < theta * max( 1, upper ) )
         break
      end
   end
end

% The pseudo-hard case

% Find eigen decomposition and the minimal eigenvalue.
[ V, D ]   = eig( H );
neigd      = neigd + 1;
[mu, imu ] = min( diag(D) );

if ( verbose )
   gamma   = abs(V(:,imu)'*g);
   disp( [' deft_funnel_solve_TR_MS : ============ pseudo hard case: gamma = ', ...
            num2str(gamma), ' ||g|| = ', num2str(norm(g))] )
end

% Compute the critical step and its orthogonal complement along the
% eigenvectors corresponding to the most negative eigenvalue
D        = D - mu * eye(n);
maxdiag  = max( diag( D ) );
ii       = find( abs(diag(D)) < 1e-10*maxdiag );
if ( length( ii ) < n & isempty( ii ) ~= 1 )
   D(ii,ii)    = 0.5 * maxdiag * eye(length(ii));
   Dinv        = inv(D);
   Dinv(ii,ii) = 0;
   scri        = -V*Dinv*V'*g;
   nscri       = norm( scri );
else
   scri  = zeros( n, 1 );
   nscri = 0;
end
if ( nscri <= Delta )
   root  = roots([norm(V(:,imu))^2 2*V(:,imu)'*scri nscri^2-Delta^2]);
   s     = scri + root(1)*V(:,imu);
else
   s     = Delta * scri / nscri;
end 
lambda   = -mu;
if ( verbose )
   disp([ ' deft_funnel_solve_TR_MS : ============ ||scri|| = ',            ...
          num2str( norm(scri) ), ' lambda = ', num2str( lambda) ])
end
hardcase = 1;

% Compute the model value and its gradient.
w     = H * s;
value = g' * s  + 0.5 * s' * w;
gplus = g + w;
norms = norm( s );

% Define information message.
if ( norms < ( 1 - eps_D ) * Delta )
   msg = [ 'interior solution ( ', int2str( nfact ), ...
          ' factorizations,  lambda = ', num2str( lambda ), ')' ];
else
   msg = [ 'boundary solution ( ', int2str( nfact ), ' factorizations, ',   ...
           int2str( neigd ), ' eigen decomposition, lambda = ',             ...
           num2str( lambda ), ' )' ];
end

if ( verbose )
   disp( ' deft_funnel_solve_TR_MS : ============ hard case exit' )
end

end % end of deft_funnel_solve_TR_MS