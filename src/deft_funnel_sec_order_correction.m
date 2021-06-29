function [sample_set, iterate, setting, evaluations, xstatus, sstatus,      ...
    Delta_f, Delta_z, Delta_f_counter, vmax, rho, pos, it_type,             ...
    poised_model] = deft_funnel_sec_order_correction(f, c, h, sample_set,   ...
    iterate, iterate_plus, setting, models, indicators, evaluations,        ...
    d, delta_f, vmax, pos, xstatus, sstatus, const, model_size, Delta,      ...
    Delta_f, Delta_z, Deltamax, Delta_f_counter, rho, succ, it_type,        ...
    poised_model)
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Desc: Calculates a second-order correction for the normal step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = iterate.xdim;
m = iterate.sdim;

nmax = min(Delta_z, setting.kappa_n * norm(iterate_plus.z_s));
lb(1:n+m)   = -nmax;
ub(1:n+m)   = nmax;
dlx         = setting.lx(iterate.indfree) - iterate_plus.x;
dlx         = dlx';
dux         = setting.ux(iterate.indfree) - iterate_plus.x;
dux         = dux';
lb(1:n)     = max(lb(1:n), dlx);
ub(1:n)     = min(ub(1:n), dux);

dls         = setting.ls - iterate_plus.s;
dls         = dls';
dus         = setting.us - iterate_plus.s;
dus         = dus';
lb(n+1:n+m) = max(lb(n+1:n+m), dls);
ub(n+1:n+m) = min(ub(n+1:n+m), dus);

full_n  = length(iterate.xfix);
I = eye(full_n);

% Compute the second-order correction
sol      = deft_funnel_blls(models.derivatives.J_s, -iterate_plus.z_s, lb', ub');
n_xsoc   = sol(1:n);
n_ssoc   = sol(n+1:n+m);
norm_soc = norm(d+sol);

% Initialize the 'second-order-correction' point
iterate_soc = iterate_plus;

if (norm_soc <= Delta)

    iterate_soc.x = iterate_plus.x + n_xsoc ;
    iterate_soc.s = iterate_plus.s + n_ssoc ;
    
    % Evaluate the functions at the trial point without including it in X yet
    [sample_set, iterate_soc, evaluations, xstatus, sstatus] =              ...
        deft_funnel_eval_functions(f, c, h, iterate_soc, sample_set,        ...
        setting, evaluations, xstatus, const.unused, sstatus, 'eval_only');
    
    iterate_soc.z_s = iterate_soc.zeval - iterate_soc.s;
    vsoc = 0.5 * (iterate_soc.z_s.' * iterate_soc.z_s);

    rho = (iterate.feval - iterate_soc.feval + setting.rho_eps) /           ...
        (delta_f + setting.rho_eps);

    if (rho >= setting.eta1 && vsoc <= vmax)

        Delta_f_counter = 0;

        % Set index of new point
        sample_set.nbPoints = sample_set.nbPoints + 1;

        % Set initial status of the new point
        xstatus(sample_set.nbPoints) = 0;
        % The point is contained in current subspace
        sstatus(sample_set.nbPoints) = 1;

        % Augment X with the new point
        if (iterate_soc.nfix > 0)
            I       = eye(iterate_soc.fulldim);
            xfix    = iterate_soc.xfix;
            indfix  = iterate_soc.indfix;
            indfree = iterate_soc.indfree;
            sample_set.X(:, sample_set.nbPoints) = I(:,indfix) *            ...
                xfix(indfix) + I(:,indfree) * iterate_soc.x;
        else
            sample_set.X(:, sample_set.nbPoints) = iterate_soc.x;
        end

        if (strcmp(setting.type_f, 'BB'))
            % Augment fX with new function value
            if (max(size(iterate_soc.feval)) > 1)
                sample_set.fX(sample_set.nbPoints) = min(setting.fxmax,     ...
                    norm(iterate_soc.feval)^2);
            else
                sample_set.fX(sample_set.nbPoints) = min(setting.fxmax,     ...
                    real(iterate_soc.feval));
            end
        end

        if (setting.cons_c)
            % Augment cX with new function value
            length_c = max(size(iterate_soc.ceval));
            if (length_c > 1)
                for i=1:length_c
                    sample_set.cX(i, sample_set.nbPoints) =                 ...
                        min(setting.fxmax, real(iterate_soc.ceval(i)));
                end
            else
                sample_set.cX(sample_set.nbPoints) = min(setting.fxmax,     ...
                    real(iterate_soc.ceval));
            end
        end

        if (setting.verbose >=2)
            disp(' Successful f-iteration with 2nd-order correction :')
            disp(' ')
        end

        it_type = 'f-succ(2nd order correction used)';
        
        [sample_set, iterate, setting, pos, xstatus, Delta_f,              	...
            Delta_z, Delta_f_counter, vmax, poised_model] =                	...
            deft_funnel_succ_iteration(1, sample_set, iterate,             	...
            iterate_soc, setting, indicators, model_size, Delta_f,         	...
            Delta_z, Delta_f_counter, Deltamax, vsoc, vmax,                 ...
            xstatus, const, succ, norm_soc, norm_soc, rho, poised_model);
    end
end
    
end % end of deft_funnel_sec_order_correction
