function [ s, lambda, norms, value, gplus, nfact, neigd, msg ] =            ...
            deft_funnel_solve_TR_MS_bc( g, H, lb, ub, Delta, eps_D, stratLam )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: An implementation of exact trust-region minimization based on the
% Moré-Sorensen algorithm subject to bound constraints.
%
% Input: 
%   - g        : model's gradient
%   - H        : model's Hessian
%   - lb       : lower bounds on the step
%   - ub       : upper bounds on the step
%   - Delta    : trust-region's radius
%   - eps_D    : accuracy required on the equation ||s|| = Delta for a
%                boundary solution
%   - stratLam : strategy to adjust lambda to find an active bound
%                (1 - Newtons method, 0 - bisection method)
%
% Output:
%   - s        : trust-region step
%   - lambda   : Lagrange multiplier corresponding to the trust-region constraint
%   - norms    : norm of s
%   - value    : value of the model at the optimal solution
%   - gplus    : value of the model's gradient at the optimal solution
%   - nfact    : number of Cholesky factorizations required
%   - neigd    : number of eifenvalue decompositions required
%   - msg      : information message
%
% Dependecies  : deft_funnel_solve_TR_MS
% Programming  : A. Troeltzsch, S. Gratton, July 2009. 
% ( This version 14 I 2010 )
%
% CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose    = 0;
theta      = 1.0e-13;          % accuracy of the interval on lambda
eps_bound  = 1.0e-5;           % the max distance | bound - x | < eps_bound for 
                               % a boundary solution
                                 
% Initialization 
msg        = '';
lambda     = 0;                % initial lambda
value      = 0;                % initial model value
nfact      = 0;                % initialize in case of error
neigd      = 0;                % initialize in case of error
Delta0     = Delta;
g0         = g;                % copy initial gradient
gplus      = g;                % initialize gplus
s          = zeros(size(g));   % initial zero step
norms      = 0;
n          = length(g);        % space dimension
I          = eye(n);           % identity matrix used for projection
ind_active = [];               % indeces of active bounds
ind_free   = 1:n;              % indeces of inactive bounds
nfree      = n;                % nbr of inactive indeces

if ( length( find( isnan(H) ) ) ~= 0 )
   if (verbose )
      disp( 'Error in deft_funnel_solve_TR_MS_bc: H contains NaNs!' )
   end
   msg = 'Error1';
   return;
end

if ( length( find( ~isreal(H) ) ) ~= 0 )
   if (verbose )
      disp( 'Error in deft_funnel_solve_TR_MS_bc: H contains imaginary parts!' )
   end
   msg = 'Error2';
   return;
end

if ( length( find( isinf(H) ) ) ~= 0 )
   if (verbose )
      disp( 'Error in deft_funnel_solve_TR_MS_bc: H contains infinite elements!' )
   end
   msg = 'Error3';
   return;
end

% Fix active components
ind_g_crit = find(( abs( lb ) <= 1e-10 & g > 0 ) | ( ub <= 1e-10 & g < 0 ));

if ( length(ind_g_crit) ~= 0)
    ind_active = ind_free(ind_g_crit);
    ind_free   = setdiff(ind_free,ind_active);
    nfree      = length(ind_free);
end

% Loop until no free variables anymore
j=0;
while nfree > 0

   % Loop until all active variables are detected and fixed
   new_call_to_MS = 1;
   
   while ( new_call_to_MS == 1 )
      j=j+1;
         
      % Minimize system in the (possibly) reduced subspace
      if ( verbose >= 1 )
         disp(['(' num2str(j) ') ---- minimizing in the (sub)space of ' num2str(length(ind_free)) ' variable(s)'])
      end
        
      g_reduced = g(ind_free);
      H_reduced = H(ind_free,ind_free);

      % Call unconstrained MS in (possibly) reduced subspace
      [ s_deltaMS, lambda, norms_deltaMS, value_red, gplus_red, nfact_r,    ...
          neigd_r, msg, hardcase ] = deft_funnel_solve_TR_MS( g_reduced,    ...
          H_reduced, Delta, eps_D );
      nfact = nfact + nfact_r;
      neigd = neigd + neigd_r;
      gplus(ind_free) = gplus(ind_free) + gplus_red;

      s_after_reduced_ms = s + I(:,ind_free) * s_deltaMS;

      % Compute critical components which became active during the last 
      % MS iteration
      ind_u_crit = find((ub(ind_free)-s_after_reduced_ms(ind_free)) <= eps_bound & ub(ind_free) <= 1e-10 );
      ind_l_crit = find((s_after_reduced_ms(ind_free)-lb(ind_free)) <= eps_bound & lb(ind_free) >= -1e-10 );
        
      % Fix these new active components
      if ( length(ind_u_crit) + length(ind_l_crit) ~= 0)
           
         ind_active = [ind_active ind_free(ind_u_crit) ind_free(ind_l_crit)];
         ind_free   = setdiff(1:n,ind_active);
         nfree      = length(ind_free);
         if (verbose)
            disp('fixed one or more variables')
         end
            
         % If no inactive variables anymore --> exit
         if ( nfree == 0 )
            norms = norm( s );
            value = 0.5 * s' * H * s + s' * g0;
            if ( verbose )
               disp('no inactive variables anymore - return')
            end
            return
         end
            
      else            
         new_call_to_MS = 0;
      end
   end
    
   % Check if step is outside the bounds
   if (verbose==2)
      disp('check if step inside bounds')
   end

   out_of_ubound = find( (ub(ind_free)-s_after_reduced_ms(ind_free)) < 0.0 );
   out_of_lbound = find( (s_after_reduced_ms(ind_free)-lb(ind_free)) < 0.0 );
   out_of_ubound_init = out_of_ubound;
   out_of_lbound_init = out_of_lbound;

   if ( length(out_of_ubound) + length(out_of_lbound) ~= 0 )

      back_inside = 0;
      lambda0 = lambda;
      if (verbose==2)
         disp( 'step outside bounds!' )
         out_of_ubound
         out_of_lbound
         disp( [ 'lambda_0=' num2str(lambda0) ] )
      end

      % Set lower bound on lambda.
      lower = lambda;
        
      % New lambda for bisection method
      if ( stratLam == 0 )
         lambda = max( 2.0, 2 * lambda );
      end
        
      % Compute upper bound on lambda (using the closest bound out of the hit bounds) 
      gnorm = norm( g );
      if ( length(out_of_ubound) > 0 )
         delta_b = min( abs(ub(ind_free(out_of_ubound))-s(ind_free(out_of_ubound))) );
      end
      if ( length(out_of_lbound) > 0 )
         delta_b = min( abs(lb(ind_free(out_of_lbound))-s(ind_free(out_of_lbound))) );
         if ( length(out_of_ubound) > 0 )
            delta_b  = min( min( abs(ub(ind_free(out_of_ubound))-s(ind_free(out_of_ubound))) ), delta_b);
         end
      end
        
      goverD   = gnorm / delta_b;
      Hnorminf = norm( H, inf );
      if ( Hnorminf > 0 )        % because Octave generates a NaN for null matrices.
         HnormF   = norm( H, 'fro' );
      else
         HnormF   = 0;
      end   
      upper    = max( 0, goverD + min( Hnorminf, HnormF ) );

      % Compute active components
      ind_u_active = find(abs(ub(ind_free)-s_after_reduced_ms(ind_free)) <= eps_bound);
      ind_l_active = find(abs(s_after_reduced_ms(ind_free)-lb(ind_free)) <= eps_bound);

      % Loop on successive trial values for lambda.
      i = 0;
      while ( ( ( length(ind_u_active) + length(ind_l_active)) == 0 ) ||    ...
              ( length(out_of_lbound) + length(out_of_ubound) ~= 0 ) )
         i = i + 1;

         % Print lambda value
         old_lambda = lambda;
         new_lambda = -1;
         if ( verbose )
            disp( [ ' deft_funnel_solve_TR_MS_bc (',int2str(i),'): lower = ',  ...
                    num2str( lower ), ' lambda = ', num2str( lambda ) ,     ...
                    ' upper = ',  num2str( upper ) ] )
         end

         % Factorize H + lambda * I.
         [ R, p ] = chol( H(ind_free,ind_free) + lambda * eye( nfree ) );

         if ( length( find( isnan( R ) ) ) ~= 0 )
            H
            lambda
            norm(g)
            R
            p
            disp( 'Error in deft_funnel_solve_TR_MS_bc: NaNs in Cholesky factorization' )
            msg = 'Error4';
            return
         end
         nfact    = nfact + 1;

         % Successful factorization 
         if ( p == 0 & hardcase == 0 )

            warning off
            s_deltaH = - R \ ( R' \ g(ind_free) );

            s_duringH = s + I(:,ind_free) * s_deltaH;

            % Find components which are at its bound and became active
            ind_u_crit  = find( (ub(ind_free)-s_duringH(ind_free)) <= eps_bound & ub(ind_free) <= 1e-10 );
            ind_l_crit  = find( (s_duringH(ind_free)-lb(ind_free)) <= eps_bound & lb(ind_free) >= -1e-10 );

            % Set these new active components to zero for one iteration
            if ( length(ind_u_crit) ~= 0)
               s_deltaH(ind_u_crit)  = 0.0;
               s_duringH(ind_free(ind_u_crit)) = 0.0;
            end
            if ( length(ind_l_crit) ~= 0)
               s_deltaH(ind_l_crit)  = 0.0;
               s_duringH(ind_free(ind_l_crit)) = 0.0;
            end

            out_of_ubound = find( (ub(ind_free)-s_duringH(ind_free)) < 0.0 );
            out_of_lbound = find( (s_duringH(ind_free)-lb(ind_free)) < 0.0 );

            % Find an appropriate bound for the next homotopy-step when using Newton's method
            if ( stratLam == 1 || verbose > 0 )

               if (( length(out_of_ubound) ~= 0 ) || ( length(out_of_lbound) ~= 0 ))

                  % OUTSIDE the bounds: find the furthest step component 
                  outside = 1;
                  if ( length(out_of_ubound) ~= 0 )
                     [diff_b_u,ind_b_u] = max(abs(ub(ind_free(out_of_ubound))-s_duringH(ind_free(out_of_ubound))));
                     norms_b = abs(s_deltaH(out_of_ubound(ind_b_u)));
                     delta_b = abs(ub(ind_free(out_of_ubound(ind_b_u)))-s(ind_free(out_of_ubound(ind_b_u))));
                     ind_b   = out_of_ubound(ind_b_u);
                     sign_b  = sign(ub(ind_free(out_of_ubound(ind_b_u)))-s(ind_free(out_of_ubound(ind_b_u))));
                     out_of_ubound_init = [ out_of_ubound_init; out_of_ubound ];
                  end
                  if ( length(out_of_lbound) ~= 0 )
                     [diff_b_l,ind_b_l] = max(abs(s_duringH(ind_free(out_of_lbound))-lb(ind_free(out_of_lbound))));
                     norms_b = abs(s_deltaH(out_of_lbound(ind_b_l)));
                     delta_b = abs(lb(ind_free(out_of_lbound(ind_b_l)))-s(ind_free(out_of_lbound(ind_b_l))));
                     ind_b   = out_of_lbound(ind_b_l);
                     sign_b  = sign(lb(ind_free(out_of_lbound(ind_b_l)))-s(ind_free(out_of_lbound(ind_b_l))));
                     out_of_lbound_init = [ out_of_lbound_init; out_of_lbound ];
                  end
                  if (( length(out_of_ubound) ~= 0 ) && ( length(out_of_lbound) ~= 0 ))
                     if (diff_b_u > diff_b_l)
                        norms_b = abs(s_deltaH(out_of_ubound(ind_b_u)));
                        delta_b = abs(ub(ind_free(out_of_ubound(ind_b_u)))-s(ind_free(out_of_ubound(ind_b_u))));
                        ind_b   = out_of_ubound(ind_b_u);
                        sign_b  = sign(ub(ind_free(out_of_ubound(ind_b_u)))-s(ind_free(out_of_ubound(ind_b_u))));
                     else
                        norms_b = abs(s_deltaH(out_of_lbound(ind_b_l)));
                        delta_b = abs(lb(ind_free(out_of_lbound(ind_b_l)))-s(ind_free(out_of_lbound(ind_b_l))));
                        ind_b   = out_of_lbound(ind_b_l);
                        sign_b  = sign(lb(ind_free(out_of_lbound(ind_b_l)))-s(ind_free(out_of_lbound(ind_b_l))));
                     end
                  end

               else
                  % INSIDE the bounds but no step component active:
                  % find the closest components to its bound from the
                  % set of components which were initially outside
                  outside = 0;
                  
                  if ( length(out_of_ubound_init) ~= 0 )
                     [diff_b_u,ind_b_u] = min(abs(ub(ind_free(out_of_ubound_init))-s_duringH(ind_free(out_of_ubound_init))));
                     norms_b = abs(s_deltaH(out_of_ubound_init(ind_b_u)));
                     delta_b = abs(ub(ind_free(out_of_ubound_init(ind_b_u)))-s(ind_free(out_of_ubound_init(ind_b_u))));
                     ind_b   = out_of_ubound_init(ind_b_u);
                     sign_b  = sign(ub(ind_free(out_of_ubound_init(ind_b_u)))-s(ind_free(out_of_ubound_init(ind_b_u))));
                  end
                  if ( length(out_of_lbound_init) ~= 0 )
                     [diff_b_l,ind_b_l] = min(abs(s_duringH(ind_free(out_of_lbound_init))-lb(ind_free(out_of_lbound_init))));
                     norms_b = abs(s_deltaH(out_of_lbound_init(ind_b_l)));
                     delta_b = abs(lb(ind_free(out_of_lbound_init(ind_b_l)))-s(ind_free(out_of_lbound_init(ind_b_l))));
                     ind_b   = out_of_lbound_init(ind_b_l);
                     sign_b  = sign(lb(ind_free(out_of_lbound_init(ind_b_l)))-s(ind_free(out_of_lbound_init(ind_b_l))));
                  end
                  if (( length(out_of_ubound_init) ~= 0 ) && ( length(out_of_lbound_init) ~= 0 ))
                     if ( diff_b_u < diff_b_l )
                        norms_b = abs(s_deltaH(out_of_ubound_init(ind_b_u)));
                        delta_b = abs(ub(ind_free(out_of_ubound_init(ind_b_u)))-s(ind_free(out_of_ubound_init(ind_b_u))));
                        ind_b   = out_of_ubound_init(ind_b_u);
                        sign_b  = sign(ub(ind_free(out_of_ubound_init(ind_b_u)))-s(ind_free(out_of_ubound_init(ind_b_u))));
                     else
                        norms_b = abs(s_deltaH(out_of_lbound_init(ind_b_l)));
                        delta_b = abs(lb(ind_free(out_of_lbound_init(ind_b_l)))-s(ind_free(out_of_lbound_init(ind_b_l))));
                        ind_b   = out_of_lbound_init(ind_b_l);
                        sign_b  = sign(lb(ind_free(out_of_lbound_init(ind_b_l)))-s(ind_free(out_of_lbound_init(ind_b_l))));
                     end
                  end
               end
            end
                 
            % Iteration printout
            if ( verbose )
               lambda_save(i)  = lambda;
               norms_b_save(i) = norms_b;

               if ( outside == 0 )
                  fprintf( 1, '%s%d%s %12.8e %s %12.8e %s\n', ' deft_funnel_solve_TR_MS_bc (', i,  ...
                  '): |s_i| = ', norms_b, '  |bound_i| = ', delta_b , '   s < bounds' )
               else
                  fprintf( 1, '%s%d%s %12.8e %s %12.8e\n', ' deft_funnel_solve_TR_MS_bc (', i, ...
                  '): |s_i| = ', norms_b, '  |bound_i| = ', delta_b )
               end
            end

            % Test if step inside bounds +/- eps_bound
            out_of_uEpsbound = find((ub(ind_free)-s_duringH(ind_free)) < -eps_bound);
            out_of_lEpsbound = find((s_duringH(ind_free)-lb(ind_free)) < -eps_bound);

            if ( isempty(out_of_uEpsbound) && isempty(out_of_lEpsbound) )
                    
               if (verbose>=2)
                  disp('all components inside the bounds + eps_bound')
               end
                    
               back_inside = 1;

               % Check if at least one component active
               ind_u_active = find(abs(ub(ind_free)-s_duringH(ind_free)) <= eps_bound);
               ind_l_active = find(abs(s_duringH(ind_free)-lb(ind_free)) <= eps_bound);

               if ( (length(ind_u_active) + length(ind_l_active)) ~= 0 )
                  if (verbose>=2)
                     disp(['all components inside the bounds + eps_bound, ' ...
                     num2str(length(ind_u_active)+length(ind_l_active))     ...
                     ' component/s close to one of its bounds'])
                  end

                  % Compute the current step after the homotopy-step
                  s_afterH = s + I(:,ind_free) * s_deltaH;

                  % Move active components to their bounds 
                  if ( length(ind_u_active) > 0 )
                     s_afterH(ind_free(ind_u_active)) = ub(ind_free(ind_u_active));
                  end
                  if ( length(ind_l_active) > 0 )
                     s_afterH(ind_free(ind_l_active)) = lb(ind_free(ind_l_active));
                  end

                  % Define information message.
                  msg = 'boundary solution';

                  break;
               end
            end

            % Compute new lambda
            if ( stratLam == 0 ) % Using bisection method
  
               if ( back_inside == 0 )
                  lambda = 2 * lambda;
                  if ( upper < lambda )
                     upper = 2*lambda; 
                  end
               else
                  if (isempty(out_of_ubound) && isempty(out_of_lbound))
                     upper = lambda;
                  else
                     lower = lambda;
                  end
                  new_lambda  = (lambda + old_lambda) / 2;
                  theta_range = theta * ( upper - lower );
                  if ( new_lambda > lower + theta_range && new_lambda < upper - theta_range )
                     lambda = new_lambda;
                  else
                     lambda = max( sqrt( lower * upper ), lower + theta_range );
                  end
               end
                     
	    else   % Using Newton's iteration on the modified secular equation

               % Reset bounds on lambda
               if (isempty(out_of_ubound) && isempty(out_of_lbound))
                  upper = lambda;
               else
                  lower = lambda;
               end
                    
               % Check the sign of the chosen step component
               unitvec        = zeros(nfree,1);
               unitvec(ind_b) = 1;
               es             = unitvec'*s_deltaH;
                    
               if ( sign(es) ~= sign_b )
                        
                  % if step component has the wrong sign
                  % (other than the active bound): one bisection step
                  new_lambda  = (lower + upper) / 2;
               else
	              % Newton step
                  w1 = R' \ unitvec;
                  w2 = R' \ s_deltaH;

                  new_lambda = lambda + ( ( norms_b - delta_b ) / delta_b ) * ( norms_b^2 / ( es*(w1'*w2) ) );
                  if ( back_inside == 0 && upper <= new_lambda )
                     upper = 2 * new_lambda;
                  end
               end

               % Check new value of lambda wrt its bounds
               theta_range = theta * ( upper - lower );
               if ( new_lambda > lower + theta_range && new_lambda <= upper - theta_range )
                  lambda = new_lambda;
               else
                  lambda = max( sqrt( lower * upper ), lower + theta_range );
               end

            end

         else % Unsuccessful factorization: compute new lambda 
                
            if ( verbose )
               disp('unsuccessful factorization')
            end
            hardcase = 0;
            lower  = lambda;
            t      = 0.5;
            lambda = ( 1 - t ) * lower + t * upper;
         end

         % Return with error message after 100 iterations
         if ( i >= 100 )
            s(1:n) = 0.0;
            norms  = 0;
            msg    = 'Error6';
            disp('Error in deft_funnel_solve_TR_MS_bc: iteration limit in bc-MS exceeded!')
            return
         end

      end % end of while-loop
            
   else
        
      % MS was successful inside the bounds
      if ( verbose >= 2 )
         disp( 'step inside bounds!' ) 
      end
        
      % Define information message.
      msg = '(partly) interior solution';
               
      % Update the step and model value
      s     = s_after_reduced_ms;
      norms = norm( s );
      value = 0.5 * s' * H * s + s' * g0;
      break
   end    
       
   % Update the step and model value
   s     = s_afterH;
   norms = norm( s );   
   value = 0.5 * s' * H * s + s' * g0;
    
   % Update trust-region radius
   Delta = Delta0 - norms;
    
   if ( Delta < -eps_bound )
      disp('Error in deft_funnel_solve_TR_MS_bc: delta smaller than zero !!!!!!')
      msg = 'Error7'; 
      return
   elseif ( Delta < 0 )
      return
   end
    
   % Update gradient 
   g = g0 + H * s;
    
   % Update active bounds
   ind_active = find( (ub-s) <= eps_bound | (s-lb) <= eps_bound );
   ind_active = ind_active';
   ind_free   = setdiff(1:n,ind_active);
   nfree      = length(ind_free);
    
   if ( nfree > 0 )
       
      % Check first-order criticality
      ng_reduced = norm( g(ind_free), inf );
      if ( ng_reduced <= 1e-5 )
         if ( verbose >= 2 )
            disp('point first order critical - return')
            ng_reduced
         end
         return
      end
    
      % Fix components which became active
      ind_g_crit = find(( abs(lb(ind_free)) <= 1e-10 & g(ind_free) > 0 ) | ( ub(ind_free) <= 1e-10 & g(ind_free) < 0 ));

      if ( length(ind_g_crit) ~= 0)
         ind_active = [ind_active ind_free(ind_g_crit)];
         ind_free   = setdiff(ind_free,ind_active);
         nfree      = length(ind_free);
      end
   end

end

end % end of deft_funnel_solve_TR_MS_bc