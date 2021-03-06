function [models, iterate, indicators] =                                    ...
    deft_funnel_update_iteration(sample_set, iterate, dev_f, dev_h, setting)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Updates the models, the indicators and the iterate structs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n  = iterate.xdim;
m  = iterate.sdim;
Im = eye(m);
In = eye(n);

% Update the models and calculate associated derivatives
models = deft_funnel_build_models(sample_set, setting);
derivatives = deft_funnel_compute_derivatives(models, sample_set, ...
    iterate, dev_f, dev_h, setting);

[mu, ws, zs, wx, zx] = deft_funnel_compute_lag_mu(iterate, derivatives, setting);

glagMatrix = [derivatives.J' zeros(n, 2 * m) In -In]; 
glagMatrix = [glagMatrix; -Im Im -Im zeros(m, 2 * n)];
LagMu = [mu; ws; zs; wx; zx];
glag = glagMatrix * LagMu;
glag = glag + [derivatives.gfx; zeros(m, 1)];

derivatives.J_s = [derivatives.J -Im];
HLag = deft_funnel_buildHLag(mu, derivatives);
derivatives.HLag = HLag;

% Add computed derivatives to models struct
models.derivatives = derivatives;

% Update the iterate struct
iterate.s = deft_funnel_slack_update(iterate.s, iterate.zeval); % does nothing for the time being
iterate.z_s = iterate.zeval - iterate.s;
iterate.mu = mu;

% Update indicators values
indicators.norm_glag = norm(glag);
[trial_pi_f, exitc_opt] = deft_funnel_compute_optimality(iterate, derivatives, setting);
if exitc_opt == 0
    indicators.pi_f = trial_pi_f;
else
    fprintf('\nWarning: Computation of the optimality measure was unsuccessful.');
    fprintf('Setting it to the norm of the Lagrangian gradient.\n');
    indicators.pi_f = indicators.norm_glag;
end
indicators.norm_z_s = norm(iterate.z_s);
indicators.v = 0.5 * (iterate.z_s.' * iterate.z_s);
indicators.pi_v = norm(derivatives.J_s' * iterate.z_s);

if (indicators.norm_z_s > 1.0e-05)
    indicators.chi_v = indicators.pi_v / indicators.norm_z_s;
else
    indicators.chi_v = 0.0;
end

% Track iterates
min_distance = 1.0e-7;
if (iterate.nfix > 0)
    
    I           = eye(iterate.fulldim);
    xfix        = iterate.xfix;
    indfix      = iterate.indfix;
    indfree     = iterate.indfree;
    iteratefull = I(:,indfix) * xfix(indfix) + I(:,indfree) * iterate.x;
   
    % Check if iterate has not changed since previous iteration
    if (norm(iterate.X(end, :) - iteratefull') > min_distance)
        iterate.X   = [iterate.X; iteratefull'];
    end
else
    if (norm(iterate.X(end, :) - iterate.x') > min_distance)
        iterate.X  = [iterate.X; iterate.x'];
    end
end

end % end of deft_funnel_update_iteration
