function [best_sol, best_feval, best_indicators, best_iterate, total_eval,   ...
    nb_local_searches, fL] = deft_funnel_multistart(f, c, h, dev_f,          ...
    dev_h, n, nb_cons_c, nb_cons_h, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative-Free Trust FUNNEL (DEFT-FUNNEL) for Global Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file implements DEFT-FUNNEL with the multistart method Multi-level 
% Single Linkage (MLSL) for global optimization. The MLSL is used to select 
% the starting points of the local searches done by the trust-region-based
% SQP algorithm called in 'deft_funnel.m'.
%
% Output: this code returns the best feasible local minimum found under a
% predefined total number of simulations.
%
% Please read the README file before embarking on this journey.
% Last update: 15/05/2019.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check the type of mandatory inputs    
if (~isnumeric(n))
    msg = [' DEFT-FUNNEL GLOBALOPT error: the number of decision',          ...
           ' variables is not numeric! Terminating.'];
    disp(msg)
    return
end

if (~isnumeric(nb_cons_c))
    msg = [' DEFT-FUNNEL GLOBALOPT error: the number of black-box',         ...
           ' constraints is not numeric! Terminating.'];
    disp(msg)
    return
end

if (~isnumeric(nb_cons_h))
    msg = [' DEFT-FUNNEL GLOBALOPT error: the number of white-box',         ...
           ' constraints is not numeric! Terminating.'];
    disp(msg)
    return
end

% Initialization
rng('shuffle', 'twister');  % seed the random nb generator based on the current time
L = [];                     % set of local minima found
fL = [];                    % set of optimal values found
S = [];                     % set of sample points
started_points = [];        % indicates which points in S have been used
fS = [];                    % objective function values of points in S
cS = [];                    % constraint values of c at the points in S
hS = [];                    % constraint values of h at the points in S
mS = [];                    % merit function values of points in S
best_sol = NaN*ones(n,1);   % best solution from all local searches
best_iterate = [];          % iterate struct of the best feasible solution
best_feval = Inf;           % best feasible objective function value
best_indicators = [];       % indicators of the best solution found
total_eval = 0;             % total number of evaluations
f_global_optimum = NaN;     % known objective function value of the global optimum
nb_global_found = 0;        % number of times the global optimum was found
nb_local_searches = 0;      % number of local searches run
maxeval = 5000*n;           % maximum number of simulations (budget)
maxeval_ls = maxeval*0.7;   % maximum number of simulations per local search
tol_feas = 1.0e-4;          % tolerance for feasibility
tol_conv = 1.0e-4;          % tolerance for declaring convergence based on the 
                            % distance to the known global optimum (not
                            % used if the global optimum is not known a
                            % priori)
type_f = 'BB';              % 'BB' if f is a black-box and 'WB' otherwise
                            % Default: 'BB'
whichmodel = 2;             % Approach to build the surrogate models:
                            % 0 = Subbasis model
                            % 1 = Frobenius-norm model
                            % 2 = minimum l2-norm (Default)
                            % 3 = regression (recommended for noisy functions)

% Set the MLSL constants
nb_of_new_samples = min(20*n, floor(0.40*maxeval)); % number of new samples per MLSL iteration
sigma = 4;
kappa = 1;                  % percentage of the best samples to consider
maxit = 20;                 % number of MLSL iterations
alpha = pi^(-0.5) * (gamma(1+n/2) * sigma)^(1/n); % MLSL formula

% Setting default bounds
lx(1:n) = -Inf;             % Default 'x' lower bounds
lx=lx';      
ux(1:n) = Inf;              % Default 'x' upper bounds
ux=ux';      
ls(1:nb_cons_c+nb_cons_h) = 0.0; % Default 's' lower bounds
ls=ls';      
us(1:nb_cons_c+nb_cons_h) = 0.0; % Default 's' upper bounds
us=us';  

% Check the parity of the variable argument list.
noptargs = size(varargin, 2);    % The number of optional arguments
if mod(noptargs, 2) > 0
   msg   = [' DEFT-FUNNEL error: the number of variable arguments',         ...
            ' must be even! Default parameters used.'];
   disp(msg);
   return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%      PROCESS ARGUMENT LIST      %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:2:noptargs

    % Lower bounds on the x variables
    if (strcmp(varargin{i}, 'lxbounds'))
        if (isnumeric(varargin{i+1}))
            if (isempty(varargin{i+1}) == 0 && ...
                 length(varargin{i+1}) == n)
                lx = varargin{i+1};
                if (size(lx, 1) == 1)
                    lx=lx';
                end
            else
                msg = [' DEFT-FUNNEL GLOBALOPT warning: lower x bounds',    ...
                       ' empty or not of length n! No bounds used.'];
                disp(msg)
            end
        else
            msg = [' DEFT-FUNNEL GLOBALOPT warning: wrong type of input',   ... 
                   ' for lower x bounds! No bounds used.'];
            disp(msg)
        end
        
    % Upper bounds on the x variables
    elseif (strcmp(varargin{i}, 'uxbounds'))
        if (isnumeric(varargin{i+1}))
            if (isempty(varargin{i+1}) == 0 && ...
                 length(varargin{i+1}) == n)
                ux = varargin{i+1};
                if (size(ux, 1) == 1)
                    ux = ux';
                end
            else
                msg = [' DEFT-FUNNEL GLOBALOPT warning: upper x bounds',    ...
                       ' empty or not of length n! No bounds used.'];
                disp(msg)
            end
        else
            msg = [' DEFT-FUNNEL GLOBALOPTOPT warning: wrong type of',      ... 
                   ' input for upper x bounds! No bounds used.'];
            disp(msg)
        end
        
    % Lower bounds on the s variables
    elseif (strcmp(varargin{i}, 'lsbounds'))
        if (isnumeric(varargin{i+1}))
            if (isempty(varargin{i+1}) == 0 && ...
                 length(varargin{i+1}) == nb_cons_c + nb_cons_h)
                ls = varargin{i+1};
                if (size(ls, 1) == 1)
                    ls = ls';
                end
            else
                msg = [' DEFT-FUNNEL GLOBALOPT warning: lower s bounds',    ...
                       ' empty or not of length "nb_cons_c + nb_cons_h"!',  ...
                       ' Setting to zero by default.'];
                disp(msg)
            end
        else
            msg = [' DEFT-FUNNEL GLOBALOPT warning: wrong type of input',   ... 
                   ' for lower s bounds! Setting to zero by default.'];
            disp(msg)
        end
        
    % Upper bounds on the s variables
    elseif (strcmp(varargin{i}, 'usbounds'))
        if (isnumeric(varargin{i+1}))
            if (isempty(varargin{i+1}) == 0 && ...
                 length(varargin{i+1}) == nb_cons_c + nb_cons_h)
                us = varargin{i+1};
                if (size(us, 1) == 1)
                    us = us';
                end
            else
                msg = [' DEFT-FUNNEL GLOBALOPT warning: upper s bounds',    ...
                       ' empty or not of length "nb_cons_c + nb_cons_h"!',  ...
                       ' Setting to zero by default.'];
                disp(msg)
            end
        else
            msg = [' DEFT-FUNNEL GLOBALOPT warning: wrong type of input',   ... 
                   ' for upper s bounds! Setting to zero by default.'];
            disp(msg)
        end
        
    % Maximum number of simulations
    elseif (strcmp(varargin{i}, 'maxeval'))
        if (isnumeric(varargin{i+1}))
            maxeval = varargin{i+1};
        else
            msg = [' DEFT-FUNNEL GLOBALOPT warning: wrong type of',         ... 
                   ' input for maximum number of simulations!',            ...
                   ' Default value will be used.'];
            disp(msg)
        end
        
    % Maximum number of simulations per local search
    elseif (strcmp(varargin{i}, 'maxeval_ls'))
        if (isnumeric(varargin{i+1}))
            maxeval_ls = varargin{i+1};
        else
            msg = [' DEFT-FUNNEL GLOBALOPT error: wrong type of input',     ... 
                   ' for maximum number of simulations per local search!',  ...
                   ' Default value will be used.'];
            disp(msg)
        end
        
    % Type of the interpolation model
    elseif (strcmp(varargin{i}, 'whichmodel'))
        if (isnumeric(varargin{i+1}))
            whichmodel = varargin{i+1};
        else
            msg = [' DEFT-FUNNEL GLOBALOPT warning: wrong type of input',   ... 
                   ' for the type of the interpolation model.',             ...
                   '  Default value (minimum l2-norm) will be used.'];
            disp(msg)
        end
        
    % Type of the objective function
    elseif (strcmp(varargin{i}, 'type_f'))
        if (ischar(varargin{i+1}))
            if (strcmp(varargin{i+1},'WB') || (strcmp(varargin{i+1},'wb')))
                type_f = 'WB';
                if (nb_cons_c == 0) % No black-box functions at all
                    msg = [' DEFT-FUNNEL error: no black-box',              ... 
                           ' function found. Please check if the',          ...
                           ' objective function is of black-box type',      ...
                           ' and if there is at least one BB constraint'];
                    disp(msg)
                    return
                end
            end
        else
            msg = [' DEFT-FUNNEL GLOBALOPT warning: wrong type of input',   ... 
                   ' for the type of the objective function.',              ...
                   ' Default value (BB) will be used.'];
            disp(msg)
        end
        
    % Ojective function value of the global optimum known by the user
    elseif (strcmp(varargin{i}, 'f_global_optimum'))
        if (isnumeric(varargin{i+1}))
            f_global_optimum = varargin{i+1};
        else
            msg = [' DEFT-FUNNEL GLOBALOPT warning: wrong type of input',   ... 
                   ' for the ojective function value of the global',        ...
                   ' optimum! Ignoring input value.'];
            disp(msg)
        end
        
    elseif (strcmp(varargin{i}, 'seed'))
        if (isnumeric(varargin{i+1}))
            seed = varargin{i+1};
            rng(seed);
        else
            msg = [' DEFT-FUNNEL GLOBALOPT warning: wrong type of input',   ... 
                   ' for the seed! Ignoring input value.'];
            disp(msg)
        end
        
    else
        msg = [' DEFT-FUNNEL GLOBALOPT warning: undefined keyword',         ...
                varargin{i}, '! Ignoring input value.'];
        disp(msg)
    end
end

% Calculate a rough approximation of the volume of the feasible region.
% If the feasible region is not bounded by a box, set the 'x' lower/upper
% bounds to -10/10 when necessary
volume_factor = 1;
for i = 1:n
    a = lx(i);
    b = ux(i);
    if isinf(a)
        a = -10;
    end
    if isinf(b)
        b = max(0,a) + 10;
    end
    volume_factor = volume_factor*(b - a);
end
alpha = alpha * volume_factor^(1.0/n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%      MULTISTART SEARCH      %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:maxit
    
    if (total_eval >= maxeval)
        disp(' TERMINATED: maximum budget reached.')
        disp(' ')
        break
    end
    
    if (total_eval + nb_of_new_samples > maxeval)
        total_eval
        maxeval
        disp(' TERMINATED: the remaining budget is not enough to continue.')
        disp(' ')
        break
    end
    
    % Sample new points and add them to S
    Snew = deft_funnel_multistart_sampling('uniform', n, nb_of_new_samples, lx, ux);
    
    nbsamples_old = size(S, 2);
    S = [S Snew];
    started_points = [started_points zeros(1, nb_of_new_samples)];
    nbsamples = size(S, 2);
    
    % Calculate their objective, constraint and merit function values
    for i = 1:size(Snew, 2)
        if (strcmp(c, 'combined'))
            try
                output = f(Snew(:, i));
                fS(nbsamples_old + i)    = output(1);
                cS(:, nbsamples_old + i) = output(2:nb_cons_c+1);
            catch
                disp(' ')
                disp(' Error: evaluation of the black box FAILED at the point');
                Snew(:, i)
            end
        else
            try
                fS(nbsamples_old + i)    = f(Snew(:, i));
            catch
                disp(' ')
                disp( ' Error: evaluation of the objective function FAILED at the point');
                Snew(:, i)
            end
            if (nb_cons_c > 0)
                try
                    cS(:, nbsamples_old + i) = c(Snew(:, i));
                catch
                    disp(' ')
                    disp(' Error: evaluation of the constraint(s) c FAILED at the point');
                    Snew(:, i)
                end
            end
        end
        if (nb_cons_h > 0)
            try
                hS(:, nbsamples_old + i) = h(Snew(:, i));
            catch
                disp(' ')
                disp( ' Error: evaluation of the constraint(s) h FAILED at the point');
                Snew(:, i)
            end
        end
    end
    total_eval = total_eval + nb_of_new_samples;
    
    if (~isempty(cS))
        blackbox_evals = cS(:, nbsamples_old+1: end);
    else
        blackbox_evals = [];
    end
    if (~isempty(hS))
        whitebox_evals = hS(:, nbsamples_old+1: end);
    else
        whitebox_evals = [];
    end
    
    [mSnew, cons_viol_Snew] = deft_funnel_multistart_merit_function(Snew,   ...
        fS(nbsamples_old+1: end), blackbox_evals, whitebox_evals,           ...
        nb_cons_c, nb_cons_h, ls, us);
    mS = [mS mSnew];
    
    [Sreduced, mSreduced, asc_order] = deft_funnel_multistart_best_samples(S, mS, kappa);
    
    % Update cluster radius
    rk = alpha * (log(k*nb_of_new_samples)/(k*nb_of_new_samples))^(1/n);
        
    for j = 1:size(Sreduced, 2)
        
        if (total_eval >= maxeval)
            disp(' TERMINATED: maximum budget reached.')
            disp(' ')
            break
        end
        
        belong_cluster = deft_funnel_multistart_cluster(Sreduced, mSreduced, j, rk);
        
        % Check if there is a better point in a nearby cluster or if it
        % has already been used as a starting point
        if (~belong_cluster && ~started_points(asc_order(j)))
            
            started_points(asc_order(j)) = 1;
            x0 = Sreduced(:,j);
            total_eval = total_eval - 1;
            nb_local_searches = nb_local_searches + 1;
            
            try 
                
                % Start local search
                [x, fx, mu, indicators, evaluations, iterate, exit_algo] =  ...
                    deft_funnel(f, c, h, dev_f, dev_h, x0,                  ...
                    nb_cons_c, nb_cons_h,                                   ...
                    'lxbounds', lx,                                         ...
                    'uxbounds' , ux,                                        ...
                    'lsbounds', ls,                                         ...
                    'usbounds' , us,                                        ...
                    'multistart_call', nb_local_searches,                   ...
                    'maxeval', min(maxeval_ls, maxeval - total_eval),       ...
                    'type_f', type_f,                                       ...
                    'whichmodel', whichmodel);
                
            catch
                disp(' ')
                disp([' Local search number ', int2str(nb_local_searches), ' FAILED.'])
                disp(' Trying the next starting point.');
                exit_algo = 1;
                evaluations.nfeval = 0;
            end
              
            if (exit_algo == 0) % no errors found
                % Add new local optimum
                L = [L x];
                fL = [fL fx];
                if (fx < best_feval && indicators.norm_z_s <= tol_feas)
                    best_feval = fx;
                    best_sol = x;
                    best_indicators = indicators;
                    best_iterate = iterate;
                elseif (~isempty(iterate.best_feas) &&                      ...
                        iterate.best_feas_feval < best_feval)
                    best_feval = iterate.best_feas_feval;
                    best_sol = iterate.best_feas;
                    best_indicators = indicators;
                    best_iterate = iterate;
                end
                
                if (~isnan(f_global_optimum))
                    if (abs(fx - f_global_optimum) < tol_conv)
                        nb_global_found = nb_global_found + 1;
                    end
                end
            end
           
            total_eval = total_eval + evaluations.nfeval;          
        end
    end
end

end % end of deft_funnel_multistart
