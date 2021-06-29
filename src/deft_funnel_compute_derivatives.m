function derivatives = deft_funnel_compute_derivatives(models, sample_set,  ...
    iterate, dev_f, dev_h, setting)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Computes the derivatives of the surrogate models.
%
% Input:
%   - models      : struct of the surrogate models
%   - sample_set  : struct of the sample set
%   - iterate     : struct of pthe current iterate
%   - dev_f       : function handle to compute the derivatives of f if any
%   - dev_h       : function handle to compute the derivatives of h if any
%   - setting     : struct of parameters
%
% Output:
%   - derivatives : struct of the derivatives at the current iterate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gradient and hessian of the objective function
if (strcmp(setting.type_f, 'BB'))
    gfx = deft_funnel_gradP(models.f, iterate.x);
    Hf  = deft_funnel_hessP(models.f, iterate.x);
else
    nfix = iterate.nfix;
    
    if (nfix > 0)
        full_n  = iterate.fulldim;
        I       = eye(full_n);
        xfix    = iterate.xfix;
        indfix  = iterate.indfix;
        indfree = iterate.indfree;
        fullsp_iterate = I(:,indfix) * xfix(indfix) + I(:,indfree) * iterate.x;
    else
        fullsp_iterate = iterate.x;
    end
    
    output_dev_f = dev_f(fullsp_iterate);
    gfx = output_dev_f{1};
    if (size(gfx,2) > 1)
        gfx = gfx';   
    end
    Hf  = output_dev_f{2};
    
    % Consider only the indices corresponding to the current subspace
    if (nfix > 0)
        gfx = gfx(indfree);
        Hf = Hf(indfree, indfree);
    end
end

% Initialize the hessian of the constraints
J = [];
HZlist = [];
    
% Compute the Jacobian and the Hessian of the contraint functions c
if (setting.cons_c)
    for i=1:setting.nb_cons_c
        J(i,:) = deft_funnel_gradP(models.c(i,:), iterate.x);
        Ci = deft_funnel_hessP(models.c(i,:), iterate.x);
        HZlist = [HZlist Ci];
    end
end
    
% Add the Jacobian and the Hessian of the contraint functions h
if (setting.cons_h)
    nfix = iterate.nfix;

    if (nfix > 0)
        full_n  = iterate.fulldim;
        I       = eye(full_n);
        xfix    = iterate.xfix;
        indfix  = iterate.indfix;
        indfree = iterate.indfree;
        fullsp_iterate = I(:,indfix) * xfix(indfix) + I(:,indfree) * iterate.x;
    else
        fullsp_iterate = iterate.x;
    end
    
    output_dev_h = dev_h(fullsp_iterate);
    Jh = output_dev_h{1};
    Hhlist = output_dev_h{2};
    
    % Consider only the indices corresponding to the current subspace
    if (nfix > 0)
        Jh = Jh(:, indfree);
        Hhlist_subsp = [];
        for i=1:setting.nb_cons_h
            Hh_i = Hhlist(:, (i-1)*full_n+1:i*full_n);
            Hh_i = Hh_i(indfree, indfree);
            Hhlist_subsp = [Hhlist_subsp Hh_i];
        end
        Hhlist = Hhlist_subsp;
    end
    
    J = [J; Jh];
    HZlist = [HZlist Hhlist];
    
end
    
% Assemble them all into a single structure
derivatives.gfx = gfx;
derivatives.Hf = Hf;
derivatives.J = J; 
derivatives.HZlist = HZlist;

end % enf of deft_funnel_compute_derivatives
