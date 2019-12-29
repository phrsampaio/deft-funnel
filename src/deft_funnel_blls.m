function [s, resn, exitc] = deft_funnel_blls(A, b, lb, ub, searchMethod)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLLS is a solver for the bound-constrained linear least-squares problem
%
%          min || As - b ||^2     subject to    lb <= s <= ub,
%
% where ||.|| is the Euclidean norm and the contraint is understood 
% componentwise.  No restriction is made on the dimension of A or its
% rank, except that n>0 and m>0.
%
% Input: 
%  - A    :  n x m matrix
%  - b    :  m x 1 vector b
%  - lb   :  n x 1 vector of (possibly -Inf) lower bounds lb
%  - ub   :  n x 1 vector of (possibly +Inf) upper bounds ub
%
% Output:
%  - s    : (approximate) solution
%  - resn : residual norm, i.e. ||As-b||
%  - exitc:  0: s contains the problem solution
%            1: s contains the best approximate solution found in 
%               maxiter = 2*n iterations
%           -1: the algorithm failed in the sense that the Cauchy step (see
%               could not be computed in maxbck backtackings (should not happen)
%           
% The method is intended for small-dimensional problems. It is an active-set
% algorithm where the unconstrained problem is solved at each iteration in
% the subspace defined by the currently active bounds, themselves being
% determined by a projected Cauchy step. Each subspace solution is computed
% using a SVD decomposition of the reduced matrix.
%
% Programming: Ph. Sampaio, Ph. Toint, A. Troeltzsch, April 2014.
% Last modified: 21/06/2019.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m, n]   = size(A);       % problem size

% Names
armijo  = 1;              % Armijo bactracking search
spwmin  = 2;              % successive piecewise minimization search

% Algorithmic parameters 
verbose = 0;              % verbosity flag (0: silent, larger -> more output)
epsfeas = 1.0e-12;        % the tolerance for bound feasibility (primal)
epsconv = 1.0e-10;        % the tolerance for optimality        (dual)
epsres  = 1.0e-12;        % The tolerance for zero residual
maxiter = 2*n;            % the maximum number of iterations
epsdzer = 1.0e-14;        % the tolerance for considering a component of
                          % a search direction to be zero
armijob = 0.5;            % the backtracking factor in the Armijo search
armijor = 0.01;           % the linear reduction factor in the Armijo search
maxback = 15;             % the maximum number of backtracks in the Armijo
                          % searches for and beyond the Cauchy point. This
                          % value and armijob should be chosen such that 
                          % armijob^maxback is of the order of the desired 
                          % precision on the variables.

if nargin < 5
    searchMethod = 1;
end

cauchy = searchMethod;
subsrch = searchMethod;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
inds  = [1:n];
nit   = 0;                % the iteration counter
nuns  = 0;                % the number of unsuccessful bactrackings (if
                          % set to 1, disables backtracking beyond the
			              % Cauchy point)

ssub  = zeros(n, 1);      % prepare the shape of the subspace solution
exitc = 0;                % initialize exit condition to successful

% Start from the projection of the unconstrained minimizer 
% on the feasible set.
s     = min(max(pinv(A) * b, lb), ub);

% Compute the associated residual and optimality measure
res   = A*s - b;                             % initial residual    
resn  = norm(res);                           % initial residual norm
g     = A'*res;                              % initial gradient
stry  = min(max(s-g, lb), ub) - s;
opt   = norm(stry);                          % initial optimality meas

atlb      = find(abs(lb - s) <= epsfeas)';   % inds of vars at lower bound
atub      = find(abs(ub - s) <= epsfeas)';   % inds of vars at upper bound
atb       = [atlb atub];
latlb     = length(atlb);   % the number of variables at their lower bound
latub     = length(atub);   % the number of variables at their upper bound
free      = inds;
free(atb) = [];
lfree     = length(free);   % the number of free variables

% Print output banner and information on the initial iterate
if (verbose > 0)
   disp(' ')
   disp('   **************************************************************')
   disp('   *                                                            *')
   disp('   *                          BLLS                              *')
   disp('   *                                                            *')
   disp('   *   a direct bound-constrained linear least-squares solver   *')
   disp('   *                                                            *')
   disp('   *                                                            *')
   disp('   *   (c) Ph. R. Sampaio, Ph. L. Toint, A. Troeltzsch, 2014    *')
   disp('   *                                                            *')
   disp('   **************************************************************')
   disp(' ')
   disp(['     The problem has ',int2str(n),' variables and ', ...
         int2str(m), ' rows.'])
   disp(' ')
   if (verbose > 2)
      problem_matrix  = A
      right_hand_side = b'
      lower_bounds    = lb'
      upper_bounds    = ub'
      disp(' ')
   end
   fprintf('   nit         ||r||     optimality')
   fprintf('               nfr nlow nupp\n\n' )
   fprintf('  %4d       %.4e  %.4e              %4d %4d %4d\n', ...
            nit, resn, opt, lfree, latlb, latub)
      
   if (verbose > 1)
      if (verbose > 2)
          unconstrained_solution = s'
      end
      disp(' ')
      disp('   --------------------------------------------------------------')
      disp(' ')
   end
end

% Terminate if the true solution is the initial one
if ((opt <= epsconv && resn <= epsres*1.0e+7) || resn <= epsres)
   maxit = 0;   % Skips useless iterations and branches to final printout
else
   maxit = maxiter;
end

% Main iteration loop
for i = 1:maxit

   nit = nit + 1;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %        Compute a Cauchy point using a projected Armijo search.         %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   g = A' * res; % the quadratic's gradient

   % Two strategies are available for the computation of the Cauchy point:
   % 1) a simple Armijo backtracking method starting from the stepsize
   %    corresponding to the unconstrained minimizer along the steepest descent;
   % 2) a successive piecewise quadratic minimization, where the quadratic
   %    model is successively minimized on each segment of the projected
   %    steepest-descent path until a (first) local minimizer is found.

   %%%%%%%%%%%%%%%%%% Use successive piecewise minimization  %%%%%%%%%%%%%%%%

   if (cauchy == spwmin)

      % Print Cauchy search banner
      if (verbose > 1)
         disp(['   Iteration ', int2str(i)])
         disp(' ')
         fprintf('   Cauchy point search')
         fprintf(' by piecewise minimization\n')
         fprintf('          i                quad        linear')
         fprintf('     nfr nlow nupp\n')
      end

      % Recompute the set of free, lower-bounded and upper-bounded variables
      atlb = atlb(find(g(atlb) > 0));
      atub = atub(find(g(atub) < 0));
      atb  = [ atlb atub ];
      free = inds;
      free(atb) = [];

      % Perform the successive minimization
      [s, res, free, natlb, natub] =                                        ...
           deft_funnel_blls_spwmin(s, -g, A, g, res, lb, ub, free, atlb,    ...
           atub, verbose);

      % Update the variables' activity status and the residual norm
      lfree = length(free);
      latlb = length(atlb);
      latub = length(atub);
      resn  = norm(res);

   %%%%%%%%%%%%%%%%%%%%%%%%% Use Armijo backtracking  %%%%%%%%%%%%%%%%%%%%%%%

   elseif (cauchy == armijo)

      % Print Cauchy search banner
      if ( verbose > 1 )
         disp(['   Iteration ', int2str(i) ] )
         disp(' ' )
         fprintf('   Cauchy point search by projected bactracking\n')
         fprintf('     k  t      ||r||      stepsize')
         fprintf('                nfr nlow nupp\n')
         fprintf('  %4d       %.4e                          %4d %4d %4d',   ...
                  0, resn, lfree, latlb, latub)
      end

      % Perform the Armijo search
      alpha = (norm(g)/norm(A*g))^2;               % the minimizing stepsize
      rtry  = res;
      for kc = 1:maxback                           % Cauchy point loop
          stry = min(max(s-alpha*g, lb), ub);
          dtry = stry - s;
          ltry = g'*dtry;
          Adir = A*dtry;
          qtry = ltry + 0.5*Adir'*Adir;            % the new quadratic value
          rtry = rtry + Adir;
          if (verbose > 1)
             fprintf('\n  %4d  p    %.4e  %.4e', kc, norm(A*stry-b), alpha)
          end 
          if (qtry <= armijor*ltry)                % the Armijo condition
             break
          else
             alpha = armijob * alpha;              % backtracking
          end
      end % of the Armijo search loop

      % Exit if the Cauchy point calculation fails (this should not happen!)
      if (kc >= maxback)
         exitc = -1;
         break
      end
      
      % Determine the bounds which are active at the Cauchy point
      s         = stry;
      res       = rtry;
      resn      = norm(res);
      atlb      = inds(find(abs(s-lb)<= epsfeas));
      atub      = inds(find(abs(s-ub)<= epsfeas));
      atb       = [atlb atub];
      latlb     = length(atlb);
      latub     = length(atub);
      free      = inds;
      free(atb) = [];
      lfree     = length(free);
   
   end
   
   if (verbose > 1)
      if (cauchy == armijo)
         fprintf('              %4d %4d %4d\n', lfree, latlb, latub)
      end
      if (verbose > 2)
         Cauchy_point                              = s'
         indices_of_free_variables                 = free
         indices_of_variables_at_their_lower_bound = atlb
         indices_of_variables_at_their_upper_bound = atub
      end
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %       Beyond the Cauchy point: nested subspace minimization            %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % The Cauchy point is optimal for the problem because either there is no 
   % free variable left (which implies opt = 0) or because the residual is 
   % zero: terminate as s is the problem's solution.
   if (lfree ==  0 || resn <= epsres)
            
      if (verbose > 1)
         fprintf('   No nested subspace search\n')
      end
      opt  = 0;       %  Forces termination

   % There are free variables at the Cauchy point: start a nested 
   % subspace search
   else  % ( lfree > 0 && resn > epsres )

      % Print initial information for nested subspace search, if requested
      if (verbose > 1)
         fprintf('   Nested subspace search')
         if ( subsrch == spwmin )
            fprintf(' by piecewise minimization\n' )
            fprintf('     k  t      ||r||                ')
         elseif (subsrch == armijo)
            fprintf(' using projected bactracking\n' )
            fprintf('     k  t      ||r||      stepsize  ')
         end
         fprintf('              nfr nlow nupp\n')
         fprintf('  %4d       %.4e                          %4d %4d %4d\n', ...
                 0, resn, lfree, latlb, latub)
      end

      % Loop on the successive nested subspaces
      for k = 1:n

         % Solve in the subspace of currently free variables
         if (verbose > 2) 
            disp(['    > Solving in subspace ', int2str(k)])
            indices_of_free_variables                 = free
            indices_of_variables_at_their_lower_bound = atlb
            indices_of_variables_at_their_upper_bound = atub
         end
     
         if (latlb > 0)
            if (latub > 0)
               ssub(free) = pinv(A(1:m,free)) * ...
                     (b - A(1:m,atlb)*lb(atlb) - A(1:m,atub)*ub(atub));
            else
               ssub(free) = pinv(A(1:m,free)) * (b - A(1:m,atlb)*lb(atlb));
            end
         elseif (latub > 0)
            ssub(free) = pinv(A(1:m,free)) * (b - A(1:m,atub)*ub(atub));
         else
            ssub(free) = pinv(A(1:m,free)) * b;
         end

         % Build the full-space version of the subspace minimizer
         % and compute the associated residual and its norm
         ssub(atlb) = lb(atlb);
         ssub(atub) = ub(atub);
         rsubo      = A * ssub - b;
         rsubon     = norm(rsubo);

         % Check feasibility of the subspace solution wrt to inactive bounds
         natlb = find(ssub(free) < lb(free));
         natub = find(ssub(free) > ub(free));
         lnatb = length(natlb) + length(natub);

         % If the subspace minimizer is unfeasible, attempt a projected
         % backtracking Armijo search until (at worst) back in the 
         % current subspace
         if (lnatb > 0)
  
            alpha  = 1;  
            dtry   = ssub   - s;        % the direction to the free minimizer

            % Two strategies are available for the minimization of the model
            % quadratic along the "direction" dtry:
            % 1) a simple Armijo bactraking method starting from the unit 
            %   stepsize, possibly followed (in case of failure) by a step to
            %   the closest constraint boundary,
            % 2) a successive piecewise quadratic minimization, where the 
            %    quadratic model is successively minimizer on each segment 
            %    of the projected path until a (first) local minimizer is found.

            %%%%%%%%%%% Use a successive piecewise minimization  %%%%%%%%%%%%

            if (subsrch == spwmin)

               if (verbose > 1)
                  fprintf('          i                quad        linear\n')
               end

               g       = zeros(n,1);
               g(free) = A(1:m,free)' * res;      % the free gradient
               [s, res, free, atlb, atub] = ...
                       deft_funnel_blls_spwmin(s, dtry, A, g, res, ...
                                    lb, ub, free, atlb, atub, verbose);

               % Update the variables' activity status and the residual norm
               lfree = length(free);
               latlb = length(atlb);
               latub = length(atub);
               resn  = norm(res);

               if (verbose > 1)
                  fprintf('  %4d  m    %.4e             ', k, resn)
                  fprintf('             %4d %4d %4d\n', lfree, latlb, latub)
                  if (verbose > 2)
                     current_subspace_solution                 = ssub'
                     indices_of_free_variables                 = free
                     indices_of_variables_at_their_lower_bound = atlb
                     indices_of_variables_at_their_upper_bound = atub
                     free
                   end
                end

            %%%% Use Armijo backtracking + step to next violated bound  %%%%%

            elseif (subsrch == armijo)

               rred   = rsubon - resn;     % maximal residual reduction
               found  = 0;                 % backtracking success indicator

               % Backtracking loop
               nback = 4*(1-nuns);         % completely heuristic here!
                                           % any better idea?

               for kb = 1:nback

                  % Compute the unprojected trial point
                  stry = (1 - alpha) * s + alpha * ssub; 

                  % See which of the inactive bound are violated 
                  % at the unprojected trial point
                  natlbt = free(find(stry(free) < lb(free)));
                  natubt = free(find(stry(free) > ub(free)));
                  lnatbt = length(natlbt)+length(natubt);

                  % Project the trial point
                  stry = min(max(stry, lb), ub);

                  % Print information on the new point, if requested
                  if (verbose > 0)
                     rtry        = A * stry - b;
                     rtryn       = norm(rtry);
                     atlbt       = [atlb natlbt];
                     atubt       = [atub natubt];
                     atbt        = [atlbt atubt];
                     freet       = inds;
                     freet(atbt) = [];
                     latlbt      = length(atlbt);
                     latubt      = length(atubt);
                     lfreet      = length(freet);
                     fprintf('  %4d  p    %.4e  %.4e', kb, rtryn, alpha)
                     fprintf('              %4d %4d %4d\n', lfreet, latlbt, ...
                         latubt)
                  end

                  % If the unprojected trial point is in the current subspace, 
                  % terminate backtracking
                  if (lnatbt == 0)
                     break
                  end

                  % If sufficient reduction in residual norm is obtained,
                  % the computed point is a suitable next iterate.  
                  % Branch out after computing the residual and the activities
                  if (verbose == 0)
                     rtry  = A * stry - b;
                     rtryn = norm(rtry);
                  end
                  if (rtryn <= resn - armijor * alpha * rred)
                     s     = stry;
                     res   = rtry;
                     resn  = rtryn;
                     if (verbose == 0)
                        atlb      = [atlb natlbt];
                        atub      = [atub natubt];
                        atb       = [atlb atub];
                        free      = inds;
                        free(atb) = [];
                        latlb     = length(atlb);
                        latub     = length(atub);
                        lfree     = length(free);
                     else
                        atlb      = atlbt;
                        atub      = atubt;
                        free      = freet;
                        latlb     = latlbt;
                        latub     = latubt;
                        lfree     = lfreet;
                     end
                     found = 1; % Current point ok as next iterate
                     break
                  end

                  % Prepare for one more step of backtracking
                  alpha = armijob * alpha;

               end % of the bactraking loop

               % A suitable iterate has been found by backtracking: no further
               % action is needed.
               if (found)
                  break

               % The backtracking search failed to return a suitable iterate:
               % compute the maximum feasible step toward the minimizer in the
               % subspace determined by the free variables at the Cauchy point.
               else

                  if (kb >= nback)
                     nuns = nuns + 1;             % Increment failure counter
                  end
                  alpha  = 1;
                  for kf = 1:length(free)         % Compute max step along dtry
                     kk = free(kf);
                     if (dtry(kk) >= epsdzer)
                        alpha = min(alpha, (ub(kk)-s(kk))/dtry(kk));
                     elseif (dtry(kk) <= - epsdzer)
                        alpha = min( alpha, (lb(kk)-s(kk))/dtry(kk));
                     end
                  end
                  ssub   = s + alpha * dtry;      % Take max step along dtry
                  rsub   = (1 - alpha) * res + alpha * rsubo; % Residual
                  rsubn  = norm(rsub); 

                  % Print information on the new point, if requested
                  if (verbose > 1)
                     fprintf('  %4d  s    %.4e  %.4e', k, rsubn, alpha)
                     fprintf('              %4d %4d %4d\n', ...
                              lfree, latlb, latub)
                  end

                  % Determine which new bounds are active and add them to 
                  % the current active set
                  natlb     = free(find(abs(ssub(free)-lb(free))<= epsfeas));
                  natub     = free(find(abs(ssub(free)-ub(free))<= epsfeas));
                  atlb      = [atlb natlb];
                  atub      = [atub natub];
                  atb       = [atlb atub];
                  free      = inds;
                  free(atb) = [];
                  latlb     = length(atlb);
                  latub     = length(atub);
                  lfree     = length(free);

                  % Print detailed activity information for the new subspace 
                  % iterate, if requested
                  if (verbose > 2)
                     current_subspace_solution                 = ssub'
                     indices_of_free_variables                 = free
                     indices_of_variables_at_their_lower_bound = atlb
                     indices_of_variables_at_their_upper_bound = atub
                  end

                  % Prepare the next subspace minimization by memorizing the
                  % current solution
                  s    = ssub;
                  res  = rsub;
                  resn = rsubn;

               end % of the bactracking + maxstep strategy

            end 

         % If the subspace minimizer is feasible, accept it as the next iterate
         else
            s    = ssub;
            res  = rsubo;
            resn = rsubon;
            if (verbose > 1)
               fprintf('  %4d  f    %.4e             ', k, resn)
               fprintf('             %4d %4d %4d\n', lfree, latlb, latub)
               if (verbose > 2)
                  current_subspace_solution                 = ssub'
                  indices_of_free_variables                 = free
                  indices_of_variables_at_their_lower_bound = atlb
                  indices_of_variables_at_their_upper_bound = atub
               end
            end
            break
         end

      end % of the loop on nested subspaces

      % Compute the optimality measure
      opt = norm(min(max(s-A'*res, lb), ub) - s);

   end % of branch on the number of free variables at the Cauchy point

   % Iteration printout
   if (verbose == 1)
      fprintf( '  %4d       %.4e  %.4e              %4d %4d %4d\n', ...
               nit, resn, opt, lfree, latlb, latub)
   elseif (verbose > 1)
      disp(' ')
      fprintf('   nit         ||r||     optimality')
      fprintf('               nfr nlow nupp\n\n')
      fprintf('  %4d       %.4e  %.4e              %4d %4d %4d\n', ...
               nit, resn, opt, lfree, latlb, latub)
      if (verbose > 2)
         current_solution                          = s'
         indices_of_free_variables                 = free
         indices_of_variables_at_their_lower_bound = atlb
         indices_of_variables_at_their_upper_bound = atub
      end
      disp('   --------------------------------------------------------------')
      disp(' ')
   end

   % Optimality termination test
   if ((opt <= epsconv && resn <= epsres*1.0e+7) || resn <= epsres)
      break
   end

end % of the main iteration loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Finish off before returning                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% See if the maximum number of iterations has been reached without other
% problem. If yes, set the exit condition accordingly.
if (exitc == 0 && nit >= maxiter)
   exitc = 1;
end

% Termination printout
   
if (verbose > 0)
   disp(' ')
   if (verbose > 2)
      indices_of_free_variables                 = free
      indices_of_variables_at_their_lower_bound = atlb
      indices_of_variables_at_their_upper_bound = atub
      final_solution                            = s'
      final_residual                            = res'
   end
   if (exitc == 1) 
      disp( '   !!! maxit reached !!!' )
   elseif (exitc == -1) 
      disp('   !!! Cauchy point calculation failure :-(  !!!' )
   else
      disp('   ---> Solved.')
   end
   disp(' ')
end

end % end of deft_funnel_blls_exp