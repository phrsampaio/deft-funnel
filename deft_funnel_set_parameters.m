function setting = deft_funnel_set_parameters( n, nbcons, cur_degree, rep_degree, initial_Y )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Defines the majority of the parameters' setting. A description
% of each parameter is given below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setting.alpha            = 0.1;            % Subsequent reduction in gradient norm 
                                           % before repair
setting.beta             = 1;              % Ratio between gradient norm and radius 
                                           % after repair
setting.CNTsin           = 0;              % Variable to produce a quasi-random vector
setting.criterion_S      = 'weighted';     % Selection of outgoing point at 
                                           % successful iterations:
setting.criterion_FP     = 'distance';     % The same, but for far points at 
                                           % unsuccessful iterations
setting.criterion_CP     = 'standard';     % The same, but for close points at 
                                           % unsuccessful iterations
setting.cur_degree       = cur_degree;     % The degree for the initial model
setting.rep_degree       = rep_degree;     % The minimum degree for the model 
                                           % after repair
setting.epsilon          = 1.0e-4;         % Gradient termination accuracy
setting.epsilon0         = 1.0e-2;         % Initial gradient accuracy threshold
setting.eps_TR           = 0.0000001;      % Rel. accuracy on the trust-region 
                                           % constraint for steps
setting.factor_CV        = 1;              % Constant for termination that 
                                           % multiplies epsilon
setting.factor_FPU       = 10;             % Multiple of TR radius defining far points
                                           % (for unsuccessful iterations)
setting.factor_FPR       = 10;             % Multiple of TR radius defining far points
                                           % (for reconstruction of a poised 
                                           % interpolation set
                                           % for large model gradient)
setting.fmax             = 1.e25;          % Maximum possible function value. Otherwise,
                                           % it is considered as NaN
setting.rho_eps          = 1.0e-14;        % Used for calculating ´rho´
setting.eps_bnd     = setting.epsilon/10;  % Threshold to define a bound as 
setting.eps_L            = 0.001;          % Rel. accuracy on the trust-region 
                                           % constraint for L max
setting.eta1             = 1.0e-4;         % Min rho value for successful iterations
setting.eta2             = 0.9;            % Min rho value for very successful iterations
setting.eta3             = 0.6 ;           % Threshold for theta used to change 
                                           % Delta_c in a f-iteration 
setting.factor_Dmax      = 1.e5;           % The ratio between Deltamax and Delta0
setting.factor_fmax      = 1.e20;          % The ratio between the upper bound 
                                           % on the objective function
                                           % value and the absolute value of t(f(x0))
setting.gamma1           = 0.01;           % Lower bound on Delta interval 
                                           % (unsuccessful its)
setting.gamma2           = 0.5 ;           % Upper bound on Delta interval 
                                           % (unsuccessful its)
setting.gamma3           = 2.0;            % Increase in Delta (successful its)
setting.gamma4           = 8.0;            % Increase in Delta (very successful its)
setting.hardcons         = 0;              % Apply hard bounds when maximizing the 
                                           % Lagrange polynomials
setting.initial_Y        = initial_Y;      % Geometry of the initial interpolation set
                                           % 'simplex' or 'random'
setting.kappa_b          = 0.9 ;           % Condition to calculate tanget step: 
                                           %      || n || <= kappa_b * Delta
setting.kappa_ca         = 1.0e+2 ;        % For the initialization of vmax
setting.kappa_cr         = 2.0 ;           % For the initialization of vmax
setting.kappa_delta 	 = 0.2 ;           % Condition to be a f-iteration: 
                                           % delta_f >= kappa_delta*delta_ft
setting.kappa_ill        = 1e+15;          % Threshold to declare a system matrix 
                                           % as ill-conditioned
setting.kappa_n          = 1.0e+2 ;        % Feasibility problem trust region: 
                                           % || n || <= min(Delta_c, kappa_n*pi_v)
setting.kappa_tx1 	     = 0.9 ;           % vmax update formula
setting.kappa_tx2	     = 0.5 ;           % vmax update formula
setting.kappa_th         = 2000;           % Threshold for a safely nondegenerate 
                                           % set of points
setting.Lambda_XN        = 1.0e-10;        % Poisedness for new iterates
setting.Lambda_CP        = 1.5;            % Poisedness for close points
setting.Lambda_FP        = 1.0e-10;        % Poisedness for far points

setting.lSolver          = 1;              % local solver for minimizing the 
                                           % Lagrange polynomial functions
setting.maxeval          = 500 * n;        % Maximum number of evaluations
setting.maxit            = 100*setting.maxeval;% Maximum number of iterations
setting.multistart_call  = 0;              % Multistart call number
setting.scalefacX        = ones(1,n);      % Scaling factors initialized to one
setting.scaleX           = 0;              % Scaling of variables is applied
setting.show_errg        = 1;              % Display of the gradient error estimate
setting.shrink_Delta     = 1;              % Shrink trust-region radius in every 
                                           % unsuccessful iteration. 
                                           % Default: 1 (true)
setting.shrink_Delta_max = 0*n;              % Max nb of times to shrink on 
                                           % unsuccesful iterations
setting.stallfact_subsp  = 1.0e-2;            % Termination (stalled) when
                                           % ||d|| or Delta <=
                                           % epsilon*stallfact*norm(x)
setting.stallfact_fullsp = 1.0e-4;           
                                           % or Delta <= epsilon*stallfact*norm(x)
setting.stratLam         = 1;              % Strategy to adjust lambda when solving 
                                           % the bc MS problem
                                           % nearly-active: |x - bnd|<eps_bnd
setting.track_error      = 1;              % Compute the error on the models
                                           % based on the Lambda-poisedness
                                           % measure
setting.verbose          = 1;              % Printout quantity
setting.whichmodel       = 2;              % Approach to build the local models:
                                           % 0 = Subbasis model
                                           % 1 = Frobenius-norm model
                                           % 2 = minimum l2-norm
                                           % 3 = regression
                                           
% Set the default bounds for variables x and s
setting.lx(1:n) = -Inf;                    % Default x lower bounds
setting.lx=setting.lx';      
setting.ux(1:n) = Inf;                     % Default x upper bounds
setting.ux=setting.ux';      
 
setting.ls(1:nbcons) = 0.0;                % Default s lower bounds
setting.ls=setting.ls';      
setting.us(1:nbcons) = 0.0;                % Default s upper bounds
setting.us=setting.us';  

end % end of deft_funnel_set_parameters