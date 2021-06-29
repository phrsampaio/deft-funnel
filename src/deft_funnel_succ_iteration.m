function [sample_set, iterate, setting, pos, xstatus, Delta_type,           ...
    Delta_z, Delta_type_counter, vmax, poised_model] =                      ...
    deft_funnel_succ_iteration(iterType, sample_set, iterate, iterate_plus, ...
    setting, indicators, model_size, Delta_type, Delta_z,                   ...
    Delta_type_counter, Deltamax, vplus, vmax, xstatus, const, succ,        ...
    norm_vector, norm_d, rho, poised_model)
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Desc: Defines the updating steps for the iterate, the trust regions radii 
% and the interpolation set at succesfull f- and c-iterations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Update the interpolation set by including the new point
[sample_set, setting, xstatus, pos] =                                       ...
    deft_funnel_update_succ_iter_Y(sample_set, iterate_plus,                ...
    setting, model_size, xstatus, const, succ);

% If xplus could/should be included in the interpolation set,
% update the information of the new iterate
if (pos > 0)

    % Move successful trial point to the first position
    sample_set = deft_funnel_swap_in_Y(1, pos, sample_set, setting);

    % Update the iterate data
    iterate = iterate_plus;
    sample_set.i_xbest = sample_set.nbPoints;

    if (iterType == 2) % z-iteration
        
        % Update vmax
        vmax = max(setting.kappa_tx1 * vmax,                                ...
                    setting.kappa_tx2 * indicators.v +                      ...
                    (1.0 - setting.kappa_tx2) * vplus);
    else
        
        % Maybe update Delta_z
        if (vplus < setting.eta3 * vmax)
            %Delta_z = min(1.5*Delta_z,Deltamax);
            Delta_z = min(max(Delta_z, setting.gamma3 * norm_d), Deltamax);
        end
    end

    % Update Delta_type
    if (rho >= setting.eta2)
        
        if (norm_vector > 0.5*Delta_type)
            Delta_type = min(max(Delta_type, ...
                                   setting.gamma4 * Delta_type), Deltamax);
        else
            Delta_type = min(max(Delta_type, ...
                                   setting.gamma3 * Delta_type), Deltamax); 
        end
    end
    
    % Set the counter of the number of reductions of Delta_type
    % at unsuccessful iters to zero
    Delta_type_counter = 0;
    
    poised_model = 0;
end
    
end % end of deft_funnel_succ_iteration
