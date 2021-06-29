function setting = deft_funnel_set_parameters(n, nb_cons_c, nb_cons_h,      ...
    cur_degree, rep_degree, initial_Y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Defines the parameters setting.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TRUST REGION

setting.eps_TR           = 1.0e-6;         % Rel. accuracy on the trust-region 
                                           % constraint for steps
setting.rho_eps          = 1.0e-14;        % Used for avoiding rounding error 
                                           % in the calculation of 'rho'
setting.eta1             = 1.0e-4;         % Min rho value for successful iterations
setting.eta2             = 0.9;            % Min rho value for very successful iterations
setting.eta3             = 0.6;            % Threshold for theta used to change 
                                           % Delta_c in a f-iteration
setting.factor_Dmax      = 1.0e5;          % The ratio between Deltamax and Delta0
setting.gamma1           = 0.01;           % Lower bound on Delta interval 
                                           % (unsuccessful its)
setting.gamma2           = 0.5;            % Upper bound on Delta interval 
                                           % (unsuccessful its)
setting.gamma3           = 2;              % Increase in Delta (successful its)
setting.gamma4           = 3;              % Increase in Delta (very successful its)
setting.shrink_Delta     = 1;              % Shrink trust-region radius in every 
                                           % unsuccessful iteration. 
                                           % Default: 1 (true)
setting.shrink_Delta_max = n;              % Max nb of times to shrink on 
                                           % unsuccesful iterations

%% CRITICALITY

setting.alpha            = 0.1;            % Reduction ratio for epislon_i
                                           % (used in the criticality step)
setting.beta             = 1;              % Ratio between optimality/feasibility
                                           % measures and the trust-region
                                           % radius (used in the
                                           % criticality step)
setting.epsilon          = 1.0e-4;         % Optimality/feasibility termination accuracy
setting.epsilon0         = 1.0e-2;         % Initial gradient accuracy threshold
setting.maxeval          = 500*n;          % Maximum number of evaluations in a local search
setting.maxit            = 100*setting.maxeval; % Maximum number of iterations in a local search
setting.stallfact_subsp  = 1.0e-4;         % Termination (stalled) when
                                           % ||d|| or Delta <= epsilon*stallfact*norm(x)
setting.stallfact_fullsp = 1.0e-4;           
                                           
%% INTERPOLATION SET AND MODELS

setting.criterion_S      = 'weighted';     % Selection criterion of outgoing point at 
                                           % successful iterations:
setting.criterion_FP     = 'distance';     % The same, but for far points at 
                                           % unsuccessful iterations
setting.criterion_CP     = 'standard';     % The same, but for close points at 
                                           % unsuccessful iterations
setting.cur_degree       = cur_degree;     % The degree for the initial model
setting.rep_degree       = rep_degree;     % The minimum degree for the model 
                                           % after reparation
setting.factor_CV        = 10;             % Constant for termination that 
                                           % multiplies epsilon
setting.factor_FPU       = 10;             % Multiple of TR radius defining far points
                                           % (for unsuccessful iterations)
setting.factor_FPR       = 10;             % Multiple of TR radius defining far points
                                           % (for reconstruction of a poised 
                                           % interpolation set
                                           % for large model gradient)
setting.hardcons         = 0;              % Apply hard bounds when maximizing the 
                                           % Lagrange polynomials. 1 if yes
                                           % and 0 otherwise. Default: 0.
setting.initial_Y        = initial_Y;      % Geometry of the initial interpolation set
                                           % 'simplex' or 'random'
setting.kappa_th         = 2000;           % Threshold for a safely nondegenerate 
                                           % set of points
setting.Lambda_XN        = 1.0e-10;        % Poisedness for new iterates
setting.Lambda_CP        = 1.5;            % Poisedness for close points
setting.Lambda_FP        = 1.0e-10;        % Poisedness for far points
setting.lSolver          = 1;              % local solver for minimizing the 
                                           % Lagrange polynomial functions
setting.stratLam         = 1;              % Strategy to adjust lambda when solving 
                                           % the bc MS problem
setting.track_error      = 1;              % Compute the error on the models
                                           % regularly based on the Lambda-poisedness measure
                                           % NOT USED. TBD
setting.whichmodel       = 2;              % Approach to build the surrogate models:
                                           % 0 = Subbasis model
                                           % 1 = Frobenius-norm model
                                           % 2 = minimum l2-norm
                                           % 3 = regression (recommended for noisy functions)
                                           
%% GENERAL CONSTANTS

setting.CNTsin           = 0;              % Variable to produce a quasi-random vector
setting.eps_bnd          = setting.epsilon;% Threshold to define a bound as nearly active
setting.eps_L            = 0.001;          % Rel. accuracy on the trust-region 
                                           % constraint for maximization of
                                           % Lagrangian function
setting.factor_fmax      = 1.e20;          % The ratio between the upper bound 
                                           % on the objective function
                                           % value and the absolute value of t(f(x0))
setting.fmax             = 1.e25;          % Maximum possible function value. Otherwise,
                                           % it is considered as NaN
setting.kappa_r          = 0.9;            % Condition to calculate tanget step: 
                                           % | n || <= kappa_r * Delta.
                                           % Default: 0.9
setting.kappa_ca         = 1.0e+2;         % For the initialization of vmax. 
                                           % Default: 1.0e+2.
setting.kappa_cr         = 1.0;            % For the initialization of vmax
setting.kappa_delta 	 = 0.2;            % Condition to be a f-iteration: 
                                           % delta_f >= kappa_delta*delta_ft
setting.kappa_ill        = 1e+15;          % Threshold to declare a system matrix 
                                           % as ill-conditioned
setting.kappa_n          = 1.0;            % Feasibility problem trust region: 
                                           % ||n_step|| <= min(Delta_c, kappa_n*pi_v)
setting.kappa_tx1        = 0.9;            % vmax update formula
setting.kappa_tx2        = 0.5;            % vmax update formula
setting.multistart_call  = 0;              % Multistart call number
setting.nb_cons_c        = nb_cons_c;      % Number of black-box constraints
setting.nb_cons_h        = nb_cons_h;      % Number of white-box constraints
setting.scalefacX        = ones(1,n);      % Scaling factors initialized to one
setting.scaleX           = 0;              % Scaling of variables is applied
setting.show_errg        = 1;              % Display of the gradient error estimate
setting.verbose          = 1;              % Printout level

if (nb_cons_c>0)
    setting.cons_c       = true;           % True if there is at least one
                                           % black-box constraint
else
    setting.cons_c       = false;
end
if (nb_cons_h>0)
    setting.cons_h       = true;           % True if there is at least one
                                           % white-box constraint
else
    setting.cons_h       = false;
end
setting.type_f           = 'BB';
                                           
%% LOWER AND UPPER BOUNDS

setting.lx(1:n) = -Inf;                    % Default x lower bounds
setting.lx=setting.lx';      
setting.ux(1:n) = Inf;                     % Default x upper bounds
setting.ux=setting.ux';      
setting.ls(1:nb_cons_c+nb_cons_h) = 0.0;   % Default s lower bounds
setting.ls=setting.ls';      
setting.us(1:nb_cons_c+nb_cons_h) = 0.0;   % Default s upper bounds
setting.us=setting.us';  

end % end of deft_funnel_set_parameters
