function sample = deft_funnel_create_sample( x, nbcons, ref_sample )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Creates a new sample structure from scratch for a given (possibly empty) 
% vector x or uses a referential sample given as input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( ~isempty(ref_sample) )
    
    % Copy only the structure of ref_sample possibly with x as main vector
    sample         = ref_sample;
    if ( ~isempty( x ) )
        sample.x   = x;
    end
    sample.s       = [];
    sample.feval   = NaN;            % Obj. function value of the current sample
    sample.ceval   = NaN;            % Constraints values of the current sample
    sample.X       = [];
    
else
    
    % Should only be used once in 'deft_funnel.m'
    n              = length( x );
    sample.x       = x;
    sample.s       = [];
    sample.xfix    = zeros( n, 1 );
    sample.nfix    = 0;              % Number of fixed values
    sample.indfix  = [];             % Indices of fixed values in the current subspace
    sample.indfree = [];             % Indices of free values in the current subspace
    sample.xdim    = n;              % Dimension of x in the current subspace
    sample.sdim    = nbcons;         % Nb of constraints
    sample.fulldim = n;              % Dimension of x in the original problem
    sample.feval   = NaN;            % Obj. function value of the current sample
    sample.ceval   = NaN;            % Constraints values of the current sample
    sample.X       = [];             % Set of 'successors' (useful for tracking
                                     % iterate points
end

end % end of deft_funnel_create_sample