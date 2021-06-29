function [x, fx] = deft_funnel_spg(x, objfun, projfun)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Implements the Nonmonotone Spectral Projected Gradient method 
% proposed by Ernesto G. Birgin, José Mario Martínez, and Marcos Raydan, in 
% "Nonmonotone Spectral Projected Gradient Methods on Convex Sets", 
% SIAM J. Optim., 10(4), 1196–1211, 2006.
%
% It minimizes objfun(x) subject to x in C, where C is a convext set.
%
% Input:
%   - x       : initial guess
%   - objfun  : objective function handle
%   - projfun : function handle that returns projection of x onto C
%
% Output:
%   - x       : best point found
%   - fx      : objective function value of the best point found
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
n              = length(x);

epsconv        = 1e-11;   % tolerance for criticality
epsprog        = 1e-12;   % tolerance used to check for lack of progress
maxit          = 5000;    % maximum number of calls to objfun
gamma          = 1e-4;    % constant in Armijo's condition
interp_type    = 0;       % type of interpolation (0: simple backtracking; 
                          % 1: quadratic; 2: cubic)
M              = 10;      % number of steps to look back in non-monotone 
                          % Armijo's condition
projection_arc = 0;       % backtracking along projection arc (default: 0)
feasible_init  = 1;       % initial point: 0 = infeasible or 1 = feasible.
bbstep         = 1;       % type of Barzilai-Borwein step (default: 1)
verbose        = 0;       % verbosity level (0: none; 1: iteration details;
                          % 2: more details).

if (verbose > 0)
   disp(' ')
   disp('  **************************************************************')
   disp('  *                Spectral Projected Gradient                 *')
   disp('  **************************************************************')
   disp(' ')
   disp(['   The problem has ',int2str(n),' variables.'])
   disp(' ')
   disp('   Iteration     Criticality     Steplength    Function value')
   disp(' ')
   disp('  --------------------------------------------------------------')
   disp(' ')
end

% Evaluate initial point
if ~feasible_init
    x = projfun(x);
end

[fx, gx]    = objfun(x);
alpha         = 1/norm(projfun(x - gx) - x, Inf) ;
Lambda        = alpha;
listFvalues   = -Inf(M, 1);

for i = 1:maxit
    
    % Criticality test
    criticality = norm(projfun(x - gx) - x, Inf);
    
    if (criticality < epsconv)
        if (verbose > 0)
            fprintf('   %10d %15.5e %15.5e %15.5e\n', i, criticality,       ...
                     Lambda, fx);
            fprintf('\n    ****');
            fprintf('  Converged: criticality condition satisfied');
            fprintf('  ****\n');
        end
    	return;
    end
    
    % Define direction
    d = -gx;
    if (~projection_arc)
        d = projfun(x + alpha*d) - x;
    end

    % Check if progress is possible
    gtd = gx'*d;
    if (gtd > -epsprog)
        if (verbose >= 1)
            fprintf('   !!! Terminated: Directional derivative is too');
            fprintf(' small !!!\n');
        end
        return;
    end
    
    % Compute trial point
    if (projection_arc)
        Lambda = alpha;
        xplus  = projfun(x + Lambda * d);
    else
        Lambda = 1;
        xplus = x + Lambda * d;
    end
    [fxplus, gxplus] = objfun(xplus);

    % Compute the maximum function value from the past 'M' iterations
    if (M == 1)
        maxFvalue = fx;
    else
        if (i <= M)
            listFvalues(i) = fx;
        else
            listFvalues = [listFvalues(2:M); fx];
        end
        maxFvalue = max(listFvalues);
    end

    % Backtracking linesearch
    totalLinesearches = 0;
    
    while (fxplus > maxFvalue + gamma * gx' * (xplus - x))
        
        totalLinesearches = totalLinesearches + 1;
        
        if (interp_type == 0)
            
            if (verbose > 1)
                fprintf('\n   Simple backtracking: halving the stepsize.\n');
            end
            Lambda = Lambda/2;
            
        elseif (interp_type == 1 && totalLinesearches < 2)
            
            if (verbose > 1)
                fprintf('\n   Quadratic interpolation\n');
            end
            Lambda = deft_funnel_quad_interp(Lambda, fx, fxplus, gtd);
        else
            if (verbose > 1)
                fprintf('\n   Cubic interpolation\n');
            end
            Lambda = deft_funnel_cub_interp(0, Lambda, fx, fxplus, gtd, gxplus'*d);
        end

        % Return if trial step is too small
        if (norm(Lambda * d, Inf) < epsprog)
            if (verbose > 0)
                fprintf('   %10d %15.5e %15.5e %15.5e\n', i, ...
                         criticality, Lambda, fx);
                fprintf('   !!! Terminated: trial step below epsprog !!!\n');
            end
            return;
        end

        % Evaluate new trial point
        if (projection_arc)
            xplus = projfun(x + Lambda * d);
        else
            xplus = x + Lambda * d;
        end
        [fxplus, gxplus] = objfun(xplus);

    end

    % Accept trial point
    fx_old = fx;
    gx_old = gx;
    x_old  = x;
    
    x  = xplus;
    fx = fxplus;
    gx = gxplus;
    
    y = gx - gx_old;
    s = x - x_old;
    
    if (bbstep == 1)
        alpha = (s'*s) / (s'*y) ;
    else
        alpha =  (s'*y) / (y'*y);
    end
    
    if (s'*y <= 0)
        alpha = 1.0e+30;
    else
        alpha = min(1.0e+30, max(1.0e-30, alpha));
    end

    if (abs(fx - fx_old) < epsprog)
        if (verbose > 0)
            fprintf('   %10d %15.5e %15.5e %15.5e\n', i, ...
                     criticality, Lambda, fx);
            fprintf('   !!! Terminated: progress too slow !!!\n');
        end
        return;
    end
    
    % Iteration printout  
    if (verbose > 0)
        fprintf('   %10d %15.5e %15.5e %15.5e\n', i, criticality, Lambda, fx);
    end
    
end

if (verbose > 0)
    fprintf('   !!! maxit reached !!!\n');
end

end % end of deft_funnel_spg