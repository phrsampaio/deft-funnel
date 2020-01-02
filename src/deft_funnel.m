function [x, fx, mu, indicators, evaluations, iterate, exit_algo] =         ...
    deft_funnel(f, c, h, dev_f, dev_h, x0, nb_cons_c, nb_cons_h, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative-Free Trust FUNNEL (DEFT-FUNNEL) for Local Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Global Optimization, check the file 'deft_funnel_multistart.m'.
% Please read the README file before embarking on this journey.
% Last update: 12/08/2019.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check the type of mandatory inputs
if (~isnumeric(x0))
    msg = [' DEFT-FUNNEL error: the starting point is not a',               ...
           ' numeric vector! Terminating.'];
    disp(msg)
    return
end

if (~isnumeric(nb_cons_c))
    msg = [' DEFT-FUNNEL error: the number of black-box',                   ...
           ' constraints is not numeric! Terminating.'];
    disp(msg)
    return
end

if (~isnumeric(nb_cons_h))
    msg = [' DEFT-FUNNEL error: the number of white-box',                   ...
           ' constraints is not numeric! Terminating.'];
    disp(msg)
    return
end

% Reset the random generator 
rng('shuffle', 'twister');

% Make sure the starting point is a column vector
if (size(x0, 1) == 1 && size(x0, 2) > 1)
    x0 = x0';
end

% Set the dimension of the problem
n = length(x0);

% Check for faulty input dimensions 
if (nb_cons_c < 0 || nb_cons_h < 0 || n < 0)
    msg = ' DEFT-FUNNEL error: The problem dimensions are faulty.';
    disp(msg)
    return
end

% Check the parity of the variable argument list.
noptargs = size(varargin, 2); % The number of optional arguments
if mod(noptargs, 2) > 0
    msg = [' DEFT-FUNNEL error: the number of variable',                    ...
           ' arguments must be even! Default parameters used.'];
    disp(msg);
    return
end

% Set some constants
const.free          = 0;               % Variable free
const.fixed         = 1;               % Variable fixed
const.alwaysfixed   = 2;               % Variable always fixed
const.in            = 1;               % Point is in the current subspace
const.out           = 0;               % Point is not in the current subspace
const.unused        = 0;               % Point is not in the current sample set Y
const.inY           = 1;               % Point is in the current sample set Y

% Define initial values
exit_algo           = 0;               % 0 if no errors were found and 1 otherwise
nit                 = 0;               % number of iterations
nitold              = 0;               % same as above but used when entering a subspace
x                   = x0;              % Current iterate
fx                  = NaN;             % Function value of current iterate
mu                  = NaN;             % Lagrange multipliers estimates
evaluations         = [];              % Struct for evaluation counters      
indicators          = [];              % Struct for optimality and feasibility measures
sample_set.X        = [];              % Set of all points already evaluated
sample_set.fX       = [];              % Objective func values of points in X
sample_set.cX       = [];              % Constraints values of points in X
sample_set.Y        = [];              % Set of sample points currently used
sample_set.fY       = [];              % Objective func values of points in Y
sample_set.cY       = [];              % Constraints values of points in Y
sample_set.nbPoints = 0;               % Total number of evaluated points, i.e. |X|.
                                       % Updated in 'deft_funnel_eval_functions' only
sample_set.ind_Y    = [];              % Indices of points of the sample set Y in X
                                       % Updated after calling 'deft_funnel_eval_functions'
sample_set.i_xbest  = NaN;             % Indice of the current iterate (i.e., Y(1)) in X
xstatus             = [];              % Indicates whether the points are in Y or not (1 or 0)
sstatus             = [];              % Indicates whether the points are in the subspace or not (1 or 0)
vstatus             = zeros(n, 1);     % Variable status (const.free or const.fixed)
sspace_save         = [];              % All subspaces already explored
xspace_save         = [];              % Indices of entering and exiting points 
                                       % corresponding to each subspace already 
                                       % explored

model_size.pquad    = ((n+1)*(n+2))/2; % Size of a fully quadratic model
model_size.pdiag    = 2*n+1;           % Size of a diagonal model
model_size.plin     = n+1;             % Size of a linear model
evaluations.nfeval  = 0;
evaluations.nceval  = 0;
evaluations.nheval  = 0;

% Set initial default trust regions radii
Delta_f0            = 1;               % Initial TR radius for the model of f
Delta_z0            = 1;               % Initial TR radius for the model of z

% Set the main parameters
setting             = deft_funnel_set_parameters(n, nb_cons_c, nb_cons_h,   ...
                          model_size.plin, model_size.plin, 'simplex'); 

% Set the 'iterate' structure
iterate = deft_funnel_create_sample(x, nb_cons_c + nb_cons_h, []);
iterate.X = iterate.x'; % For tracking iterates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%      PROCESS ARGUMENT LIST      %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:2:noptargs

    % Lower bounds on the x variables
    if (strcmp(varargin{i}, 'lxbounds'))
        if (isnumeric(varargin{i+1}))
            if (isempty(varargin{i+1}) == 0 &&                              ...
                length(varargin{i+1}) == n)
                setting.lx = varargin{i+1};
                if (size(setting.lx, 1) == 1)
                    setting.lx=setting.lx';
                end
            else
                msg = [' DEFT-FUNNEL warning: lower x bounds',              ...
                       ' empty or not of length n! No bounds used.'];
                disp(msg)
            end
        else
            msg = [' DEFT-FUNNEL warning: wrong type of input',             ... 
                   ' for lower x bounds! No bounds used.'];
            disp(msg)
        end
        
    % Upper bounds on the x variables
    elseif (strcmp(varargin{i}, 'uxbounds'))
        if (isnumeric(varargin{i+1}))
            if (isempty(varargin{i+1}) == 0 &&                              ...
                 length(varargin{i+1}) == n)
                setting.ux = varargin{i+1};
                if (size(setting.ux, 1) == 1)
                    setting.ux = setting.ux';
                end
            else
                msg = [' DEFT-FUNNEL warning: upper x bounds',              ...
                       ' empty or not of length n! No bounds used.'];
                disp(msg)
            end
        else
            msg = [' DEFT-FUNNEL warning: wrong type of input',             ... 
                   ' for upper x bounds! No bounds used.'];
            disp(msg)
        end
        
    % Lower bounds on the s variables
    elseif (strcmp(varargin{i}, 'lsbounds'))
        if (isnumeric(varargin{i+1}))
            if (isempty(varargin{i+1}) == 0 &&                              ...
                 length(varargin{i+1}) == nb_cons_c + nb_cons_h)
                setting.ls = varargin{i+1};
                if (size(setting.ls, 1) == 1)
                    setting.ls = setting.ls';
                end
            else
                msg = [' DEFT-FUNNEL warning: lower s bounds',              ...
                       ' empty or not of length "nb_cons_c + nb_cons_h"!',  ...
                       ' Setting to zero by default.'];
                disp(msg)
            end
        else
            msg = [' DEFT-FUNNEL warning: wrong type of input',             ... 
                   ' for lower s bounds! Setting to zero by default.'];
            disp(msg)
        end
        
    % Upper bounds on the s variables
    elseif (strcmp(varargin{i}, 'usbounds'))
        if (isnumeric(varargin{i+1}))
            if (isempty(varargin{i+1}) == 0 &&                              ...
                 length(varargin{i+1}) == nb_cons_c + nb_cons_h)
                setting.us = varargin{i+1};
                if (size(setting.us, 1) == 1)
                    setting.us = setting.us';
                end
            else
                msg = [' DEFT-FUNNEL warning: upper s bounds',              ...
                       ' empty or not of length "nb_cons_c + nb_cons_h"!',  ...
                       ' Setting to zero by default.'];
                disp(msg)
            end
        else
            msg = [' DEFT-FUNNEL warning: wrong type of input',             ... 
                   ' for upper s bounds! Set to zero by default.'];
            disp(msg)
        end
        
    % Check if multistart method is being used
    elseif (strcmp(varargin{i}, 'multistart_call'))
        if (isnumeric(varargin{i+1}))
            setting.multistart_call = varargin{i+1};
        else
            msg = [' DEFT-FUNNEL warning: wrong type of input',             ... 
                   ' for "multistart_call" entry!'];
            disp(msg)
        end
    
    % Maximum number of simulations
    elseif (strcmp(varargin{i}, 'maxeval'))
        if (isnumeric(varargin{i+1}))
            setting.maxeval = varargin{i+1};
        else
            msg = [' DEFT-FUNNEL warning: wrong type of input',             ... 
                   ' for maximum number of simulations. Default value',     ...
                   ' will be used.'];
            disp(msg)
        end
        
    % Type of the interpolation model
    elseif (strcmp(varargin{i}, 'whichmodel'))
        if (isnumeric(varargin{i+1}))
            setting.whichmodel = varargin{i+1};
        else
            msg = [' DEFT-FUNNEL warning: wrong type of input',             ... 
                   ' for the type of the interpolation model.',             ...
                   '  Default value (minimum l2-norm) will be used.'];
            disp(msg)
        end
        
    % Type of the objective function
    elseif (strcmp(varargin{i}, 'type_f'))
        if (ischar(varargin{i+1}))
            if (strcmp(varargin{i+1},'WB') || (strcmp(varargin{i+1},'wb')))
                setting.type_f = 'WB';
                if (setting.nb_cons_c == 0) % No black-box functions at all
                    msg = [' DEFT-FUNNEL error: no black-box',              ... 
                           ' function found. Please check if the',          ...
                           ' objective function is of black-box type',      ...
                           ' and if there is at least one BB constraint'];
                    disp(msg)
                    exit_algo = 1;
                    return
                end
            end
        else
            msg = [' DEFT-FUNNEL warning: wrong type of input',             ... 
                   ' for the type of the objective function.',              ...
                   ' Default value (BB) will be used.'];
            disp(msg)
        end
        
    else
        msg = [' DEFT-FUNNEL warning: undefined keyword',                   ...
                varargin{i}, '! Ignoring input value.'];
        disp(msg)
    end
end

% Compute the maximal TR radius
Delta_f = Delta_f0;
Delta_z = Delta_z0;
Delta = min(Delta_f, Delta_z);
Deltamax = setting.factor_Dmax * Delta;

if (setting.multistart_call <= 1)
    deft_funnel_print_banner(setting);
else
   disp(' ')
   disp([' Local search number: ', int2str(setting.multistart_call)])
end

% Check the bounds and correct intial Delta if there is no sufficient space 
% between the bounds. Shift x0 if it is out of the box.
disp(' ')
for j=1:n
    
    % Check lower and upper bounds
    if setting.lx(j) > setting.ux(j)
        msg  = [' DEFT-FUNNEL warning: Lower bound of component ',          ...
                int2str(j), ' exceeds upper bound !!'];
        disp(msg)
        exit_algo = 1;
        return
    end
    
    % Check difference between bounds
    temp(j) = setting.ux(j) - setting.lx(j);
    
    if (temp(j) < Delta + Delta)
        if (temp(j) == 0.0)
            iterate.nfix    = iterate.nfix + 1;
            iterate.indfix  = [iterate.indfix j];
            iterate.xfix(j) = setting.lx(j);
            vstatus(j)      = const.alwaysfixed;
        else
            Delta = 0.5 * temp(j);
            [Delta_f, Delta_z] = deal(Delta);
            
            disp([' Diff. between lower and upper bound of component ',     ...
                   int2str(j), ' is less than 2*Delta0! New Delta0',        ...
                   ' set to ', num2str(Delta)]);
        end
    end
    
    % Move the starting point inside the bounds if necessary
    if (iterate.x(j) < setting.lx(j))
        iterate.x(j) = setting.lx(j);
    elseif (iterate.x(j) > setting.ux(j))
        iterate.x(j) = setting.ux(j);
    end
    
end

% Scale initial point and bounds, if user-defined
if (setting.scaleX)
    for i = 1:n
        if (setting.scalefacX(i) > 0)
            iterate.x(i)  = iterate.x(i) * setting.scalefacX(i);
            setting.lx(i) = setting.lx(i) * setting.scalefacX(i);
            setting.ux(i) = setting.ux(i) * setting.scalefacX(i);
        else
            setting.scalefacX(i) = 1;
        end
    end
end

% Reset constants if some variables are fixed
if (iterate.nfix > 0)
    nfree = iterate.fulldim - iterate.nfix;
    if (nfree <= 0)
        msg = ' DEFT-FUNNEL error: No free variables. Please, enlarge search space.';
        disp(msg);
        exit_algo = 1;
        return
    end
    iterate.indfree = setdiff(1:iterate.fulldim, iterate.indfix);
    iterate.x       = iterate.x(iterate.indfree);
    iterate.xdim    = nfree;
    if (setting.cur_degree == model_size.plin)
        setting.cur_degree = nfree + 1;
    elseif (setting.cur_degree == model_size.pdiag)
        setting.cur_degree = 2*nfree+1;
    elseif (setting.cur_degree == model_size.pquad)
        setting.cur_degree = ((nfree + 1) * (nfree + 2))/2;
    end
    if (setting.rep_degree == model_size.plin)
        setting.rep_degree = nfree+1;
    elseif (setting.rep_degree == model_size.pdiag)
        setting.rep_degree = 2*nfree+1;
    elseif (setting.rep_degree == model_size.pquad)
        setting.rep_degree = ((nfree + 1) * (nfree + 2))/2;
    end
    model_size.plin  = nfree + 1;
    model_size.pdiag = 2 * nfree + 1;
    model_size.pquad = ((nfree + 1) * (nfree + 2))/2;
else
    iterate.indfree = 1:iterate.fulldim;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%      BUILDING THE INTIAL POISED SAMPLE SET      %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sample_set, iterate, evaluations, xstatus, sstatus, poised_model, msg] =   ...
    deft_funnel_build_initial_sample_set(f, c, h, sample_set, iterate,      ...
    setting, evaluations, xstatus, sstatus, const, model_size, Delta);
    
if (strcmp(msg(2:6), 'Error'))
    disp(' ')
    disp(' DEFT-FUNNEL error: Construction of initial sample set unsuccessful.')
    x = iterate.x;
    fx = iterate.feval;
    exit_algo = 1;
    return
end

% Compute the maximal objective function value
setting.fxmax = min(setting.fmax, setting.factor_fmax*max(abs(iterate.feval), 1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%      INITIALIZATION OF THE ALGORITHM      %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set initial values for the slacks
iterate.s = zeros(iterate.sdim, 1);
MIN_THRESHOLD = -1.0e+10;
MAX_THRESHOLD = 1.0e+10;
for i=1:(nb_cons_c + nb_cons_h)
    if (isfinite(setting.ls(i)) && isfinite(setting.us(i)))
        if (setting.ls(i) > MIN_THRESHOLD && setting.us(i) < MAX_THRESHOLD)
            if (iterate.zeval(i) > setting.ls(i) && iterate.zeval(i) < setting.us(i))
                iterate.s(i) = iterate.zeval(i);
            else
                iterate.s(i) = (setting.ls(i) + setting.us(i))/2;
            end
        elseif (setting.ls(i) > MIN_THRESHOLD)
            if (iterate.zeval(i) > setting.ls(i) && iterate.zeval(i) < setting.us(i))
                iterate.s(i) = iterate.zeval(i);
            else
                iterate.s(i) = setting.ls(i);
            end
        else
            if (iterate.zeval(i) > setting.ls(i) && iterate.zeval(i) < setting.us(i))
                iterate.s(i) = iterate.zeval(i);
            else
                iterate.s(i) = setting.us(i);
            end
        end
    elseif (isfinite(setting.ls(i)))
        if (setting.ls(i) > MIN_THRESHOLD)
            if (iterate.zeval(i) > setting.ls(i) && iterate.zeval(i) < setting.us(i))
                iterate.s(i) = iterate.zeval(i);
            else
                iterate.s(i) = setting.ls(i);
            end
        else
            if (iterate.zeval(i) > setting.ls(i) && iterate.zeval(i) < setting.us(i))
                iterate.s(i) = iterate.zeval(i);
            else
                iterate.s(i) = setting.us(i);
            end
        end
    else
        if (iterate.zeval(i) > setting.ls(i) && iterate.zeval(i) < setting.us(i))
            iterate.s(i) = iterate.zeval(i);
        else
            iterate.s(i) = setting.us(i);
        end
    end
end

if (norm(iterate.zeval - iterate.s) > 1.0e+3)
    iterate.s = iterate.zeval;
end

infeasls = find(iterate.zeval < setting.ls);
if (infeasls > 0)
    iterate.s(infeasls) = setting.ls(infeasls);
end

infeasus = find(iterate.zeval > setting.us);
if (infeasus > 0)
    iterate.s(infeasus) = setting.us(infeasus);
end

% Compute Lagrange multipliers estimates and the gradient of the 
% Lagrangian function at the iterate. Build the models and get
% optimality and feasibility indicators.
[models, iterate, indicators] = ...
   deft_funnel_update_iteration(sample_set, iterate, dev_f, dev_h, setting);

% Define the initial funnel radius
vmax = max(setting.kappa_ca, setting.kappa_cr * indicators.v);

% Define the initial convergence threshold
setting.epsilon_i = setting.epsilon0;

% Print the information obtained after computing the first model
deft_funnel_print_head_info(setting);
deft_funnel_printout(nit, evaluations, iterate, setting,                    ...
    indicators, vmax, 0.0, 0.0, Delta_f, Delta_z, 0.0, '', '', sample_set);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%      CALL MAIN ROUTINE      %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[exit_algo, nit, sample_set, iterate, indicators, evaluations] =            ...
    deft_funnel_main(f, c, h, dev_f, dev_h, nit, nitold, sample_set,        ...
    iterate, setting, indicators, evaluations, models, xstatus, sstatus,    ...
    vstatus, sspace_save, xspace_save, const, model_size, Delta, Delta_f,   ...
    Delta_z, Deltamax, vmax, poised_model, 'toplevel');

x = iterate.x;
fx = iterate.feval;
mu = iterate.mu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %%%%%%%%%%%%%%%%%%      END OF DEFT-FUNNEL      %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
