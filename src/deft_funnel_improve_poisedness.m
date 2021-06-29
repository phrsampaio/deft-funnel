function [sample_set, evaluations, Delta_f, Delta_z, xstatus, sstatus,      ...
    it_type, poised_model] = deft_funnel_improve_poisedness(f, c, h,        ...
    sample_set, iterate, setting, evaluations, indicators, const, xstatus,  ...
    sstatus, exit_nstep, Delta_f, Delta_z, it_type, poised_model)
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Improve the poisedness of the interpolation set and possibly update
% the trust-regions radii.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If the error in the models is not close to zero, improve the
% models by reducing the radius of the interpolation set.
if (sample_set.errg > setting.epsilon*1.0e-3 && strcmp(it_type, 'mu-iteration'))
    
    old_errg = sample_set.errg;
    old_radius = sample_set.Y_radius;
    
    radius = 0.3*sample_set.Y_radius;
    [sample_set, replaced] = deft_funnel_repair_Y(sample_set, iterate,      ...
        setting, radius);
     
    poised_model = 1;
    sample = deft_funnel_create_sample([], [], iterate);
  
    % Compute the corresponding function values.
    for i = 1:length(replaced)
    
        j = replaced(i);

        % Set index of new point and update status of the old point
        xstatus(sample_set.ind_Y(j)) = const.unused;
        sample_set.ind_Y(j) = sample_set.nbPoints+1;
        sample.x = sample_set.Y(:, j);
        
        % Update X and evaluate function at ynew
        [sample_set, sample, evaluations, xstatus, sstatus] =               ...
            deft_funnel_eval_functions(f, c, h, sample, sample_set,         ...
            setting, evaluations, xstatus, const.inY, sstatus);
        
        % Update Y function values
        if (strcmp(setting.type_f, 'BB'))
            sample_set.fY(j) = sample_set.fX(sample_set.nbPoints);
        end
        if (setting.cons_c)
            sample_set.cY(:, j) = sample_set.cX(:, sample_set.nbPoints);
        end

    end
    
    if (setting.verbose >= 2)
        if (old_radius > sample_set.Y_radius)
            disp(' Y radius reduced to improve poisedness. ')
            disp(' ')
        end
    end
    
    if (sample_set.errg < old_errg)
        
        % If poisedness was improved, allow the trust region to increase in the
        % same order of the improvement
        it_type = strcat(it_type, ': improving poisedness');
        Delta_f = Delta_f * (1 + (old_errg - sample_set.errg));
        Delta_z = Delta_z * (1 + (old_errg - sample_set.errg));
        
    elseif (sample_set.errg > setting.epsilon_i)
        Delta_f = Delta_f * 3*setting.gamma3;
        Delta_z = Delta_z * 3*setting.gamma3;
    end

elseif (strcmp(it_type, 'mu-iteration') && ...
        indicators.norm_z_s > setting.epsilon && ...
        exit_nstep == 1)
    
    % Pathological case:
    % The current point is infeasible, but BLLs was not able to
    % compute a good normal step although the error in the interpolation
    % model is small. Try reducing the trust region radius.
    Delta_z = setting.gamma2 * Delta_z;

end
    
end % end of deft_funnel_mu_iteration
