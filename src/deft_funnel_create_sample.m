function sample = deft_funnel_create_sample(x, nb_cons, ref_sample)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Creates a new sample structure from scratch for a given (possibly empty) 
% vector x or uses a referential sample given as input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isempty(ref_sample))
    
    % Copy only the structure of ref_sample possibly with x as main vector
    sample         = ref_sample;
    if (~isempty(x))
        sample.x   = x;
    end
    sample.feval   = NaN;            % Obj. function value of the current sample
    sample.ceval   = NaN;            % Constraints values of the current sample
    
else % Should only be used once in 'deft_funnel.m'
    
    % Current iterate info
    n              = length(x);
    sample.x       = x;              % x values of current iterate
    sample.s       = [];             % s values of current iterate
    sample.xfix    = zeros(n, 1);    % fixed x values
    sample.nfix    = 0;              % Number of fixed values
    sample.indfix  = [];             % Indices of fixed values in the current subspace
    sample.indfree = [];             % Indices of free values in the current subspace
    sample.xdim    = n;              % Dimension of x in the current subspace
    sample.sdim    = nb_cons;        % Nb of constraints
    sample.fulldim = n;              % Dimension of x in the original problem
    sample.feval   = NaN;            % Obj. function value of the current sample
    sample.ceval   = NaN;            % Constraints values of c at the current sample
    sample.heval   = NaN;            % Constraints values of h at the current sample
    
    % Keep record of important info
    sample.X       = [];             % Set of 'successors' (useful for tracking
                                     % iterate points
    sample.best_feas = [];           % Feasible point (defined by setting.epsilon)
                                     % with the best f value so far
    sample.best_feas_feval = Inf;    % f value of 'sample.best_feas'

end % end of deft_funnel_create_sample
