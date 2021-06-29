function [norm_tstep, tstep, tstep_x, tstep_s, d, d_x, d_s, iterate,        ...
    indicators] = deft_funnel_tangent_step(iterate, setting, indicators,    ...
    models, Delta, nstep, nstep_x, nstep_s, M, g_n)
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Computes a tangent step if there is enough 'space' left by the normal
% step and if the optimality measure is relatively large. 
% It obtains the tangent step by solving
% (1) a linear problem with 'linprog' if H is nearly a null matrix;
% (2) a nonlinear problem with the Spectral Projected Gradient method otherwise 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = iterate.xdim;
m = iterate.sdim;
norm_n = norm(nstep, Inf);

iter_nstep = iterate;
iter_nstep.x = iter_nstep.x + nstep_x;
iter_nstep.s = iter_nstep.s + nstep_s;

g_nsmall = models.derivatives.gfx + models.derivatives.HLag * nstep_x ;
    
iter_nstep_derivatives = models.derivatives;
iter_nstep_derivatives.gfx = g_nsmall;

indicators.pi_f = deft_funnel_compute_optimality(iter_nstep,                ...
    iter_nstep_derivatives, setting);
if (indicators.pi_f < 0)
    fprintf('pi_f negative. Setting it to the gradient of the Lagrangian of previous iteration.\n');
    indicators.pi_f = indicators.norm_glag;
end

if(norm_n <= setting.kappa_r * Delta)

    iterate.mu = deft_funnel_compute_lag_mu(iter_nstep,                     ...
        iter_nstep_derivatives, setting);

    % If the dual optimality measure is sufficiently large,
    % compute a suitable tangent step t.
    if (indicators.pi_f > deft_funnel_forcing(3, indicators.norm_z_s))

        Delta_within = Delta - norm_n;
        b = zeros(m, 1);

        lb(1:n+m) = -Delta_within;
        ub(1:n+m) = Delta_within;

        dlx = setting.lx(iterate.indfree) - (iter_nstep.x);
        dlx = dlx';
        dux = setting.ux(iterate.indfree) - (iter_nstep.x);
        dux = dux';
        lb(1:n) = max(lb(1:n), dlx);
        ub(1:n) = min(ub(1:n), dux);

        dls = setting.ls - (iter_nstep.s);
        dls = dls';
        dus = setting.us - (iter_nstep.s);
        dus = dus';
        lb(n+1:n+m) = min(max(lb(n+1:n+m), dls), dus);
        ub(n+1:n+m) = max(min(ub(n+1:n+m), dus), dls);

        init_point = zeros(n+m, 1);

        % If the Hessian is nearly a null matrix, use lingprog 
        % to solve the linear problem.
        if(size(find(abs(models.derivatives.HLag) < 1.0e-10), 1) == n^2)

            options = optimset('Display','off');

            [tstep, fevallinprog, exitflag] =                               ...
                linprog(g_n, [], [], models.derivatives.J_s, b, lb', ub',   ...
                init_point, options);

            tstep_x = tstep(1:n);
            tstep_s = tstep(n+1:n+m);

        else % Otherwise solve the problem using SPG.

            objfun = @(x)deft_funnel_tangent_prob(x, M, g_n);
            projfun = @(x)deft_funnel_projection(x, models.derivatives.J_s, b, lb', ub');

            tstep = deft_funnel_spg(init_point, objfun, projfun);
            tstep = min(max(lb', tstep), ub');
            tstep_x = tstep(1:n);
            tstep_s = tstep(n+1:n+m);

        end

        norm_tstep = norm(tstep, Inf);
        d = nstep + tstep;
        d_x = d(1:n);
        d_s = d(n+1:n+m);

    else
        norm_tstep = 0.0;
        d = nstep;
        d_x = d(1:n);
        d_s = d(n+1:n+m);
        tstep = zeros(n + m, 1);
        tstep_x = tstep(1:n);
        tstep_s = tstep(n+1:n+m);
    end

else
    iterate.mu = zeros(m, 1);
    norm_tstep = 0.0;
    d = nstep;
    d_x = d(1:n);
    d_s = d(n+1:n+m);
    tstep = zeros(n + m, 1);
    tstep_x = tstep(1:n);
    tstep_s = tstep(n+1:n+m);
end 

if (setting.verbose >= 2)
    disp(' Tangent step:')
    tstep_x
    tstep_s
    disp(' Full step:')
    d_x
    d_s
end
    
end % end of deft_funnel_tangent_step
