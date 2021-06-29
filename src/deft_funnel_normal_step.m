function [nstep, nstep_x, nstep_s, exit_nstep] =                            ...
    deft_funnel_normal_step(iterate, setting, indicators, models, Delta_z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Computes a normal step if infeasibility is relatively large by using 
% the BLLS solver.
%
% Dependencies: deft_funnel_blls.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if (indicators.chi_v <= setting.epsilon && indicators.norm_z_s > setting.epsilon)

    fprintf('  ***************************************************\n');
    fprintf('  *           Infeasible stationary point           *\n');
    fprintf('  ***************************************************\n');
    exit_nstep = 2;
    return;
end

n = iterate.xdim;
m = iterate.sdim;

%if ( indicators.v > 1.0e-03*indicators.pi_f )    
if (indicators.norm_z_s > setting.epsilon)
    
    nmax = min(Delta_z, setting.kappa_n * indicators.norm_z_s) ;
    lb(1:n+m) = -nmax;
    ub(1:n+m) = nmax;

    dlx = setting.lx(iterate.indfree) - iterate.x;
    dlx = dlx';
    dux = setting.ux(iterate.indfree) - iterate.x;
    dux = dux';
    lb(1:n) = max(lb(1:n), dlx);
    ub(1:n) = min(ub(1:n), dux);

    dls = setting.ls - iterate.s;
    dls = dls';
    dus = setting.us - iterate.s;
    dus = dus';
    lb(n+1:n+m) = max(lb(n+1:n+m), dls);
    ub(n+1:n+m) = min(ub(n+1:n+m), dus);
    
    [nstep, resn, exit_nstep] = ...
        deft_funnel_blls(models.derivatives.J_s, -iterate.z_s, lb', ub', 1);
    
    percentage_error = 100*(resn/indicators.norm_z_s);
    
    % If the residual is large, try another aproach
    if (percentage_error > 40 && exit_nstep == 1)

        [nstepTRIAL, resnTRIAL, exit_nstepTRIAL] = ...
            deft_funnel_blls(models.derivatives.J_s, -iterate.z_s, lb', ub', 2);
        if (resnTRIAL < resn)
            nstep = nstepTRIAL;
            resn = resnTRIAL;
            exit_nstep = exit_nstepTRIAL;
        end
    end
    
    nstep_x = nstep(1:n);
    nstep_s = nstep(n+1:n+m);
    
    if (setting.verbose >= 2)
        disp(' Normal step:')
        Delta_z
        nstep_x
        nstep_s
        resn
        exit_nstep
    end

else
    nstep = zeros(n+m, 1);
    nstep_x = nstep(1:n);
    nstep_s = nstep(n+1:n+m);
    exit_nstep = 0;
    
    if (setting.verbose >= 2)
        disp(' Normal step not computed:')
    end
end
    
end % end of deft_funnel_normal_step
