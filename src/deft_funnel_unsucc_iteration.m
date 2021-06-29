function [sample_set, iterate, setting, vmax, xstatus, it_type,             ...
    Delta_type, Delta_type_counter, Delta_z, poised_model, changedPoint,    ...
    newPoint, distNewPoint] = deft_funnel_unsucc_iteration(pos, iterType,   ...
    sample_set, iterate, iterate_plus, setting, indicators, model_size,     ...
    const, xstatus, Delta, Delta_type, Delta_z, Delta_type_counter,         ...
    Delta_max_counter, Deltamax, delta_vsq, delta_vsqn, vplus, vmax, rho,   ...
    succ, norm_vector, norm_n, it_type, poised_model)
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Desc: Defines the updating steps for the iterate, the trust regions radii 
% and the interpolation set at unsuccesfull f- and c-iterations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    
changedPoint = [];
newPoint     = [];
distNewPoint = NaN;

% The model is not fully quadratic yet: add (if possible)
% the new point to the interpolation set and recompute the model
if ((setting.cur_degree < model_size.pquad || ...
     (setting.whichmodel == 3 && setting.cur_degree < 2*model_size.pquad)) && ...
     (rho < setting.eta1 || delta_vsq < 0.2*delta_vsqn))
 
    if (setting.verbose >= 2)
        disp(' Try to augment the size of Y by adding unsuccessful point ')
        disp(' ')
    end

    [sample_set_augmented, pY] = deft_funnel_augment_Y(sample_set, iterate_plus.x, setting);
    
    %errg_augmented_Y = deft_funnel_compute_error(sample_set_augmented, iterate, setting);
    
    sample_set = sample_set_augmented;
    setting.cur_degree = pY;
    pos = pY;

    if (setting.verbose >= 2)
        disp([' Including interpolation point ',                            ...
            int2str(setting.cur_degree), ' (augm)'])
        disp(' ')
    end

    % Update status and position of the new point
    xstatus(sample_set.nbPoints) = const.inY;
    sample_set.ind_Y(pos) = sample_set.nbPoints;
    if (strcmp(setting.type_f, 'BB'))
        sample_set.fY(pos) = iterate_plus.feval;
    end
    if (setting.cons_c)
        sample_set.cY(:, pos) = iterate_plus.ceval;
    end
    
    it_type  = strcat(it_type, ': augmenting Y');

    % Shrink trust region in unsuccessful iteration
    if(setting.shrink_Delta == 1)
       if (norm_vector ~= 0.0)
           Delta_type = setting.gamma2 * norm_vector;
           Delta_type_counter = Delta_type_counter + 1;
       else
           Delta_type = setting.gamma2 * Delta_type;
           Delta_type_counter = Delta_type_counter + 1;
       end
    end

end

% Enter if the model is already fully quadratic *or*
% iterate_plus could not yet be included in the set.
% The decision to iterate_plus here depends on possibly
% eliminating another point.
if (setting.cur_degree >= model_size.pquad || pos == 0)

    if (pos == 0 && (poised_model == 0 || Delta <= setting.epsilon_i))

        % Compute the distance of the interpolation points to
        % the current iterate. (Distinguish between the badly
        % conditioned successful and the unsuccessful case!)
        dist = zeros(1, setting.cur_degree);

        if ((rho >= setting.eta1 && iterType == 1) ||                       ...
            (rho >= setting.eta1 && iterType == 2 &&                        ...
            delta_vsq >= 0.2 * delta_vsqn && norm_vector ~= 0))

            for j=1:setting.cur_degree
                if (setting.lSolver == 1)
                    dist(j) = norm(sample_set.Y(:,j) - iterate_plus.x, Inf);
                else
                    dist(j) = norm(sample_set.Y(:,j) - iterate_plus.x, Inf);
                end
            end
        else
            for j=2:setting.cur_degree
                if (setting.lSolver == 1)
                    dist(j) = norm(sample_set.Y(:,j) - iterate_plus.x, Inf);
                else
                    dist(j) = norm(sample_set.Y(:,j) - iterate_plus.x, Inf);
                end
            end
        end

        % Compute the basic distance used to define far/close points.
        % Since the calculation of the normal step makes use of
        % Delta_z rather than Delta, the trial point may be out
        % of the trust region B(x_k,Delta), possibly causing
        % a loop on unsuccessful iterations while replacing far
        % points at c-iterations. Use Delta_z to avoid this situation.
        if (iterType == 1)
            Delta_far = Delta;
        else
            Delta_far = Delta_z;
        end

        FPlength = setting.factor_FPU * (1 + setting.eps_TR) * Delta_far;

        % Replace a far interpolation point.
        if ((rho >= setting.eta1 && iterType == 1) ||                       ...
            (rho >= setting.eta1 && iterType == 2 &&                        ...
            delta_vsq >= 0.2*delta_vsqn && norm_vector ~= 0))

            % Use weighted measure, not furthest point
            criterion_FPn = 'weighted';
        else
            criterion_FPn = setting.criterion_FP;
        end

        
        [sample_set, pos] = deft_funnel_include_in_Y(sample_set,            ...
            iterate_plus.x, setting, find(dist > FPlength),                 ...
            setting.Lambda_FP, criterion_FPn, succ);

        % If a far point has been replaced by iterate_plus
        if (pos > 0)

            it_type  = strcat(it_type, ': replacing far point');

            changedPoint = pos;
            newPoint     = sample_set.Y(:, pos);
            distNewPoint = norm(sample_set.Y(:, pos) - iterate.x, Inf);

            if (setting.verbose >= 2)
                disp([' replacing interpolation point ', int2str(pos), ' (far)'])
                disp(' ')
            end

            % Update status and position of the new point
            xstatus(sample_set.ind_Y(pos)) = const.unused;
            xstatus(sample_set.nbPoints) = const.inY;
            sample_set.ind_Y(pos) = sample_set.nbPoints;
            if (strcmp(setting.type_f, 'BB'))
                sample_set.fY(pos) = iterate_plus.feval;
            end
            if (setting.cons_c)
                sample_set.cY(:, pos) = iterate_plus.ceval;
            end
            
            % Swap points if included a successful point
            if ((rho >= setting.eta1 && iterType == 1) ||                  ...
                (rho >= setting.eta1 && iterType == 2 &&                   ...
                delta_vsq >= 0.2*delta_vsqn && norm_vector ~= 0))
              
                % Move successful trial point to the first position
                sample_set = deft_funnel_swap_in_Y(1, pos, sample_set, setting);

                % Update the iterate data
                iterate = iterate_plus;
                sample_set.i_xbest = sample_set.nbPoints;

                if (iterType == 2)

                    % Update vmax
                    vmax = max(setting.kappa_tx1 * vmax, ...
                        setting.kappa_tx2 * indicators.v + (1.0 - setting.kappa_tx2) * ...
                        vplus);
                end

                if (setting.verbose >= 2)
                    disp(' swapped point to position 1')
                    disp(' ')
                end

                it_type  = strcat(it_type, '. Added point was succ.');

                % Update Delta_type
                if (rho >= setting.eta2)

                    Delta_type = min(max(Delta_type,                        ...
                        setting.gamma3 * norm_vector), Deltamax);
                end

                if (iterType == 2)

                    % Maybe update Delta_z
                    if (vplus < setting.eta3 * vmax)
                        Delta_z = min(max(Delta_z, setting.gamma3*norm_n),  ...
                        Deltamax);
                    end
                end

            else

                % Shrink trust region in unsuccessful iteration
                if(setting.shrink_Delta == 1 &&                             ...
                    Delta_type_counter < Delta_max_counter &&               ...
                    Delta_type <= setting.epsilon_i)

                    if (norm_vector ~= 0.0)
                        Delta_type = setting.gamma2 * norm_vector;
                        Delta_type_counter = Delta_type_counter + 1;
                    else
                        Delta_type = setting.gamma2 * Delta_type;
                        Delta_type_counter = Delta_type_counter + 1;
                    end
                end

            end

        end % end of 'replacing a far point'

        % Replace a close interpolation point.
        if (pos == 0)

            if ((rho >= setting.eta1 && iterType == 1) ||                   ...
                (rho >= setting.eta1 && iterType == 2 &&                    ...
                delta_vsq >= 0.2*delta_vsqn && norm_vector ~= 0))

                criterion_CPn = 'standard'; % find best improvement
            else
                criterion_CPn = setting.criterion_CP;
            end
            if ((rho >= setting.eta1 && iterType == 1) ||                	...
                (rho >= setting.eta1 && iterType == 2 &&                   	...
                delta_vsq >= 0.2 * delta_vsqn && norm_vector ~= 0))

                % try hard to include a successful point
                Lambda_CPn = 1e-15; 
            else
                Lambda_CPn = setting.Lambda_CP;
                dist(1) = 2 * FPlength; % excludes the current iterate
            end

            [sample_set, pos] = deft_funnel_include_in_Y(sample_set,        ...
                iterate_plus.x, setting, find(dist <= FPlength),            ...
                Lambda_CPn, criterion_CPn, succ);

            if (pos > 0)

                it_type  = strcat(it_type, ': replacing close point');

                changedPoint = pos;

                newPoint = sample_set.Y(:, pos);
                distNewPoint = norm(sample_set.Y(:, pos) - iterate.x);

                if (setting.verbose >= 2)
                    disp([' replacing interpolation point ', int2str(pos), ' (close)'])
                    disp(' ')
                end

                % Update status and position of the new point
                xstatus(sample_set.ind_Y(pos)) = const.unused;
                xstatus(sample_set.nbPoints) = const.inY;
                sample_set.ind_Y(pos) = sample_set.nbPoints;
                if (strcmp(setting.type_f, 'BB'))
                    sample_set.fY(pos) = iterate_plus.feval;
                end
                if (setting.cons_c)
                    sample_set.cY(:, pos) = iterate_plus.ceval;
                end
                
                % Swap points if included a successful point
                if ((rho >= setting.eta1 && iterType == 1) ||               ...
                    (rho >= setting.eta1 && iterType == 2 &&                ...
                    delta_vsq >= 0.2*delta_vsqn && norm_vector ~= 0))

                    % Move successful trial point to the first position
                    sample_set = deft_funnel_swap_in_Y(1, pos, sample_set, setting);

                    % Update the iterate data
                    iterate = iterate_plus;
                    sample_set.i_xbest = sample_set.nbPoints;

                    if (iterType == 2)

                        % Update vmax
                        vmax = max(setting.kappa_tx1 * vmax,                ...
                            setting.kappa_tx2 * indicators.v +              ...
                            (1.0 - setting.kappa_tx2) * vplus);
                    end

                    if (setting.verbose >= 2)
                        disp(' swapped point to position 1')
                        disp(' ')
                    end

                    it_type  = strcat(it_type, '. Added point was succ.');

                    % Update Delta_type
                    if (rho >= setting.eta2)

                        Delta_type = min(max(Delta_type,                  	...
                            setting.gamma3 * norm_vector), Deltamax);
                    end

                    if (iterType == 2)

                        % Maybe update Delta_z
                        if (vplus < setting.eta3 * vmax)
                            Delta_z = min(max(Delta_z,                      ...
                            setting.gamma3 * norm_n), Deltamax);
                        end
                    end

                else
                    
                    % Shrink trust region in unsuccessful iteration
                    if(setting.shrink_Delta == 1 &&                         ...
                        Delta_type_counter < Delta_max_counter &&           ...
                        Delta_type <= setting.epsilon_i)

                        if (norm_vector ~= 0.0)
                            Delta_type = setting.gamma2 * norm_vector;
                            Delta_type_counter = Delta_type_counter + 1;
                        else
                            Delta_type = setting.gamma2 * Delta_type;
                            Delta_type_counter = Delta_type_counter + 1;
                        end
                    end
                end
            end % of 'if (pos > 0)'
        end % replace a close point
    end % fully quadratic model 

    % Decrease the radius if iterate_plus hasn't been included in Y
    if (pos == 0)

        if (setting.verbose >= 2)
            disp(' decreasing the TR radius')
            disp(' ')
        end

        % Set status of the new point
        xstatus(sample_set.nbPoints) = const.unused;
        
        % Compute new trust-region radius
        if (norm_vector ~= 0.0)
            Delta_type = setting.gamma2 * norm_vector;            
        else
            Delta_type = setting.gamma2 * Delta_type;
        end
        
        it_type  = strcat(it_type, ': reducing Delta');
    end

end

if (iterType == 2)
    Delta_z = Delta_type;
end

end % end of deft_funnel_unsucc_iteration
