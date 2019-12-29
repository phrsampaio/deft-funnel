function [exit_algo, nit, sample_set, iterate, indicators, evaluations,     ...
    models, xstatus, sstatus, vstatus, sspace_save,  xspace_save,           ...
    const, model_size, Delta, Delta_f, Delta_z, vmax, msg, xlist,           ...
    trlist] = deft_funnel_main(f, c, h, dev_f, dev_h, nit, nitold,          ...
    sample_set, iterate, setting, indicators, evaluations, models,          ...
    xstatus, sstatus, vstatus, sspace_save, xspace_save, const,             ...
    model_size, Delta, Delta_f, Delta_z, Deltamax, vmax, poised_model,      ...
    level)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Main routine of the algorithm DEFT-FUNNEL. This function is called by 
% 'deft_funnel.m' and 'deft_funnel_subspace_min.m'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
msg               = '';
n                 = iterate.xdim;
nb_cons           = iterate.sdim;

norm_d            = 100;             % Initial large value to avoid entering the
                                     % Criticality step at the beginning
norm_d_x          = 0;
norm_d_s          = 0;
rho               = 0;

Delta_f_counter   = 0;
Delta_z_counter   = 0;
Delta_max_counter = setting.shrink_Delta_max;

mu_iter           = 0;               % Max number of consecutive mu-iterations allowed
radius            = setting.epsilon; % Used in the criticality step for repairing Y
it_type           = '';
succ              = 0;               % succ = 1 if rho >= setting.eta1 and 0 otherwise

nstep             = [];              % Normal step
tstep             = [];              % Tangent step

% Keep record of (x, Delta)
xlist  = [];
trlist = [];

zeros_nm = zeros(n, nb_cons);
zeros_mn = zeros(nb_cons, n);
zeros_mm = zeros(nb_cons, nb_cons);
zeros_m1 = zeros(nb_cons, 1);

exit_algo = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%      MAIN ITERATION      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = nitold+1:setting.maxit
    
    nit = nit + 1;

    % Update the best feasbile point obtained so far
    
    if (iterate.feval < iterate.best_feas_feval &&                          ...
        indicators.norm_z_s <= setting.epsilon)
        
        full_n  = iterate.fulldim;
        I       = eye(full_n);
        xfix    = iterate.xfix;
        nfix    = iterate.nfix;
        indfix  = iterate.indfix;
        indfree = iterate.indfree;
        if (nfix > 0)
            my_sample = I(:,indfix) * xfix(indfix) + I(:,indfree) * iterate.x;
        else
            my_sample = iterate.x;
        end
        
        iterate.best_feas = my_sample;
        iterate.best_feas_feval = iterate.feval;
    end

    % Counter for the loop in the criticality test
    if(strcmp(it_type, 'criticality test'))
        criticalLoop = criticalLoop + 1;
    else
        criticalLoop = 0;
    end
    
    Delta = min(Delta_f, Delta_z);
    
    xlist  = [xlist iterate.x];
    trlist = [trlist Delta];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%      SUBSPACE MINIMIZATION STEP      %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ((indicators.norm_z_s > setting.epsilon || indicators.norm_glag > setting.epsilon) && ...
         Delta  > setting.epsilon * setting.stallfact_fullsp * norm(iterate.x) && ...
         norm_d > setting.epsilon * setting.stallfact_fullsp * norm(iterate.x))

    
        [nit, sample_set, iterate, setting, indicators, evaluations,        ...
            models, xstatus, sstatus, vstatus, sspace_save, xspace_save,    ...
            Delta, Delta_f, Delta_z, vmax, msg] =                           ...
            deft_funnel_subspace_min(f, c, h, dev_f, dev_h, nit,            ...
            sample_set, iterate, setting, indicators, evaluations, models,  ...
            xstatus, sstatus, vstatus, sspace_save, xspace_save, const,     ...
            model_size, Delta, Delta_f, Delta_z, Deltamax, vmax,            ...
            poised_model, level);
        
        if (~isempty(msg))
            if (strcmp(msg(1:5), 'Error'))
                if (strcmp(level,'toplevel'))
                    disp(msg)
                end
                msg2 = 'Error found in the subspace minimization step. Stopping the run!';
                if (setting.verbose >= 1)
                    disp(msg2);
                end
                return
            end
        end
        
        % Test for the maximum number of function evaluations
        if (evaluations.nfeval >= setting.maxeval)

            if (setting.verbose >= 1)
                disp([' !!!!!!!!!!!!!!!!!!!!!!!!!!!!! MAX. EVALUATIONS',    ...
                      ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
            end
            return
        end
        
        % Test for the maximum number of iterations.
        if (k+1 >= setting.maxit)

            if (setting.verbose >= 1)
                disp([' !!!!!!!!!!!!!!!!!!!!!!!!!!!!! MAX. ITERATIONS',     ...
                      ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
            end
            return
        end
        
        it_type = 'idenActBnd';
        
        % Check if current subspace has already been explored
        if (~isempty(xspace_save))
            if (find(xspace_save(3, :) > 0))
                
                % If the past iteration was unsuccessful, the iterate
                % remains the same. Therefore no need to repair the
                % interpolation set again
                if (strcmp(level, 'toplevel') && succ)
                    if (sample_set.errg > setting.epsilon)
                        
                        xspace_save(3, :) = 0;
                        
                        % If current subspace has already been explored,
                        % repair interpolation set in smaller radius
                        maxd = 0;
                        for i = 2:setting.cur_degree
                            maxd = max(norm(sample_set.Y(:, i)-sample_set.Y(:, 1)), maxd);
                        end
                        radius = maxd * 0.5;
                        
                        if (setting.verbose > 1)
                            fprintf(' repair set and reduce TR because ');
                            fprintf(' subspace already explored\n');
                            fprintf(' new repair radius: %d', radius);
                        end
                        
                        % Check if TR has not become unreasonably small
                        if (Delta < setting.epsilon * setting.stallfact_fullsp * norm(iterate.x))
                            disp(' STOPPED: Delta became too small')
                            return
                        end
                        
                        [sample_set, replaced] =                            ...
                            deft_funnel_repair_Y(sample_set, iterate,       ...
                            setting, radius);
                        
                        poised_model = 1;
                        
                        % Compute the corresponding function values
                        
                        % Create a sample structure for evaluating the new
                        % points
                        sample = deft_funnel_create_sample([], [], iterate);
                            
                        for i = 1:length(replaced)
                            
                            j = replaced(i);
                            
                            % Set the index of the new point and update the
                            % the status of the replaced point
                            xstatus(sample_set.ind_Y(j)) = const.unused;
                            sample_set.ind_Y(j) = sample_set.nbPoints + 1;
                            
                            % Update X and evaluate the functions at Y(:, j)
                            sample.x = sample_set.Y(:, j);
                            
                            [sample_set, ~, evaluations, xstatus, sstatus] = ...
                                deft_funnel_eval_functions(f, c, h, sample, ...
                                sample_set, setting, evaluations, xstatus,  ...
                                const.inY, sstatus);
                            
                            if (strcmp(setting.type_f, 'BB' ))
                                sample_set.fY(j)    = sample_set.fX(sample_set.nbPoints);
                            end
                            if (setting.cons_c)
                                sample_set.cY(:, j) = sample_set.cX(:, sample_set.nbPoints);
                            end
                            
                            if (strcmp(msg(1:5), 'Error'))
                                if (strcmp(level,'toplevel'))
                                    disp(msg); 
                                end
                                msg3 = 'Error found in the subspace minimization step. Stopping the run!';
                                if (setting.verbose >= 1)
                                    disp(msg3);
                                end
                                return
                            end
                        end
                        
                        % Compute associated interpolation models and
                        % indicators
                        [models, iterate, indicators] =                     ...
                            deft_funnel_update_iteration(sample_set,        ...
                            iterate, dev_f, dev_h, setting);
                        
                        % Iteration summary
                        msg = ': subspace already explored. Y repaired.';
                        it_type = strcat(it_type, msg);
                        
                        % Compute maximum gradient error estimation
                        
                        sample_set = deft_funnel_poisedness_Y(sample_set, iterate, setting);
                        
                        deft_funnel_printout(nit, evaluations, iterate,     ...
                            setting, indicators, vmax, 0, 0, Delta_f,       ...
                            Delta_z, 0, it_type, '', sample_set);
                        
                        % Set new Delta
                        [Delta_z, Delta_f, Delta] = deal(radius);
                    
                    end
                    
                else
                    if (~strcmp(level, 'toplevel'))
                        % If in the subspace, return to the toplevel to repair
                        return
                    end
                end
            end
        end
    end % end of checking for active bounds

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%      CRITICALITY STEP FOR SUBLEVELS      %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Return to toplevel if convergence was attained in the sublevel 
    % (no serious check of gradient accuracy in the sublevel)
    if  (deft_funnel_criticality_check(iterate, indicators, setting,        ...
              it_type, Delta, norm_d, level) == 2)
        it_type = 'convergence in subspace';
        msg = '';
        deft_funnel_printout(nit, evaluations, iterate, setting,            ...
            indicators, vmax, norm_d_x, norm_d_s, Delta_f, Delta_z, rho,    ...
            it_type, msg, sample_set);
        return
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%      POISEDNESS IMPROVEMENT LOOP      %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check if it has been trying to improve poisedness through
    % mu-iterations for a while
    if ((strcmp(it_type, 'mu-iteration')                                 || ...
         strcmp(it_type, 'mu-iteration: improving poisedness'))          && ...
         mu_iter >= 8)
      
      it_type = 'stopped';
      
      msg = [' **********************************',                         ...
             ' Impossible to find a descent direction',                     ...
             ' ***************************************' ];
      
      deft_funnel_printout(nit, evaluations, iterate, setting,              ...
          indicators, vmax, norm_d_x, norm_d_s, Delta_f, Delta_z, rho,      ...
          it_type, msg, sample_set);

      return;
      
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%      CRITICALITY STEP FOR TOPLEVEL      %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check if optimality conditions are satisfied or if the
    % trust-region/stepsize is too small
    if  (deft_funnel_criticality_check(iterate, indicators, setting,        ...
              it_type, Delta, norm_d, level))
      
      % If the models are at their requested degree for convergence, 
      % either we have converged or the interpolation set should be 
      % repaired for a smaller epsilon_i.
      
      % Small trust region
      if (~strcmp(it_type, '') && ~strcmp(it_type, 'mu-iteration') &&       ...
           Delta <= setting.epsilon * setting.stallfact_fullsp * norm(iterate.x))
          
          criticality = 1;
          
      % Small stepsize
      elseif (~strcmp(it_type, '')                                       && ...
              ~strcmp(it_type, 'mu-iteration')                           && ...
              ~strcmp(it_type, 'criticality test')                       && ...
              norm_d <= setting.epsilon * setting.stallfact_fullsp * norm(iterate.x))
          
          criticality = 2;
      else
          criticality = 0;
      end
      
      augment = setting.rep_degree - setting.cur_degree;
      
      if (setting.verbose >= 3)
          msg = [' cur_degree = ', num2str(setting.cur_degree),             ...
                 ' rep_degree = ', num2str(setting.rep_degree) ];
          disp( msg );
      end
      
      % No additional interpolation point is required for 
      % non-terminal repair.
      if (augment <= 0)
          
          % Compute maximum gradient error estimation
          sample_set = deft_funnel_poisedness_Y(sample_set, iterate, setting);
          
          % Terminate if the solution has been found
          if (indicators.norm_z_s <= setting.epsilon && indicators.norm_glag <= setting.epsilon)
              
              if (sample_set.errg <= setting.epsilon && strcmp(level,'toplevel'))
                  
                  it_type = 'convergence';
                  msg = [' **********************************',             ...
                         ' Converged',                                      ...
                         ' ***************************************' ];
                  deft_funnel_printout(nit, evaluations, iterate, setting,  ...
                      indicators, vmax, norm_d_x, norm_d_s, Delta_f,        ...
                      Delta_z, rho, it_type, msg, sample_set);
                  return;
              end
          end
          
          % It may be possible that the criticality step was entered 
          % because optimality and feasibility were attained with
          % models that are of poor quality, in which case one has
          % criticality = 0
          if (criticality ~= 0)
              
              if (sample_set.errg <= setting.factor_CV * setting.epsilon)
                  
                  it_type = 'convergence';
                  if (criticality == 1)
                      msg = [' **********************************',         ...
                             ' Trust-region radius small',                  ...
                             ' ***************************************'] ;
                  else
                      msg = [' **********************************',         ...
                             ' Steplength small',                           ...
                             ' ***************************************'] ;
                  end
                  deft_funnel_printout(nit, evaluations, iterate, setting,  ...
                      indicators, vmax, norm_d_x, norm_d_s, Delta_f,        ...
                      Delta_z, rho, it_type, msg, sample_set);
                  return;
              end
          end
      end
      
      % Not at a solution: improve the interpolation set
      it_type = 'criticality test';
      
      % Reduce epsilon_i if repair degree is reached
      if (augment <= 0)
          setting.epsilon_i = max(max(setting.alpha * indicators.norm_z_s,  ...
              setting.alpha * indicators.pi_f), setting.epsilon);
      end
      
      % Add new interpolation points to reach the desired 
      % degree (only entered if rep_degree is higher than linear!)
      if (augment > 0)
          
          % If gradient is small, find a new point in the epsilon-environment,
          % not in Delta (distinguish between infty-norm and 2-norm 
          % local solver)
          if (indicators.norm_z_s <= setting.epsilon && indicators.norm_glag <= setting.epsilon)
              if (setting.lSolver == 2)
                  Delta = setting.epsilon/sqrt(n);
                  setting.epsilon_i = setting.epsilon/sqrt(n);
              else
                  Delta = setting.epsilon;
                  setting.epsilon_i = setting.epsilon;
              end
          end

          % Pick a random interpolation point.
          ynew = -Delta * ones(n, 1) + 2*Delta*rand(n, 1);
          [sample_set, pY] = deft_funnel_augment_Y(sample_set, ynew, setting);
          setting.cur_degree = pY;
          sample_set.ind_Y(setting.cur_degree) = pY;

          % Optimally replace it.
         if (setting.hardcons)
            ynew = deft_funnel_find_new_yj_bc(sample_set, iterate, pY, setting, Delta);
         else
            ynew = deft_funnel_find_new_yj(sample_set, pY, setting, Delta);
         end

         sample_set = deft_funnel_replace_in_Y(sample_set, ynew, pY, setting);
         replaced = [setting.cur_degree];

      % The current interpolation set has the requested degree.
      else % augment <= 0
          
         % If indicators are good, repair in epsilon-radius, else in Delta
         % (distinguish between infty-norm and 2-norm local solver)
         if (indicators.norm_z_s <= setting.epsilon                      && ...
             indicators.norm_glag <= setting.epsilon)

            if (setting.lSolver == 2)
                radius = min(0.5*setting.epsilon/sqrt(n), 0.5*radius/sqrt(n));
            else
                radius = min(0.5*setting.epsilon, 0.5*radius);
            end
         else
             radius = max(min(Delta, setting.epsilon_i), setting.epsilon);
         end
         
         [sample_set, replaced] = deft_funnel_repair_Y(sample_set, iterate, ...
             setting, radius);
          
      end
      
      poised_model = 1;
      old_errg = sample_set.errg;
      
      % Compute new gradient error estimation
      sample_set = deft_funnel_poisedness_Y(sample_set, iterate, setting);
      
      % For some reason, it cannot improve poisedness anymore
      if (old_errg <= sample_set.errg && criticalLoop > 5)
          
          if criticalLoop == 6
              [Delta_z, Delta_f, Delta] = deal(min(indicators.norm_glag, indicators.norm_z_s));
              continue
          end
          
          it_type = 'stopped';
          
          msg = [' **********************************',                     ...
                 ' Impossible to improve poisedness',                       ...
                 ' ***************************************' ];
          
          deft_funnel_printout(nit, evaluations, iterate, setting,          ...
              indicators, vmax, 0, 0, Delta_f, Delta_z, 0, it_type, msg,    ...
              sample_set);
          return;
      end
      
      % Compute the corresponding function values
      sample = deft_funnel_create_sample([], [], iterate);
       
      for i = 1:length(replaced)
          
          j = replaced(i);
          
          % Set index of new point and update status of the old point
          xstatus(sample_set.ind_Y(j)) = const.unused;
          sample_set.ind_Y(j) = sample_set.nbPoints+1;
          sample.x = sample_set.Y(:, j);
          
          % Update X and evaluate function at ynew
          [sample_set, sample, evaluations, xstatus, sstatus] =             ...
              deft_funnel_eval_functions(f, c, h, sample, sample_set,       ...
              setting, evaluations, xstatus, const.inY, sstatus);
          
          % Update Y function values
          if (strcmp(setting.type_f, 'BB' ))
            sample_set.fY(j) = sample_set.fX(sample_set.nbPoints);
          end
          
          if (setting.cons_c)
            sample_set.cY(:, j) = sample_set.cX(:, sample_set.nbPoints);
          end
          
      end
      
      % Recompute models and indicators
      [models, iterate, indicators] = ...
          deft_funnel_update_iteration(sample_set, iterate, dev_f, dev_h, setting);
      
      % Reset the trust regions radii to a multiple of the measures
      % of optimality and feasibility
      if (augment <= 0)
          Delta_f = min(setting.beta * indicators.norm_glag, Delta);
          Delta_z = min(setting.beta * indicators.norm_z_s, Delta);
          [Delta_z, Delta_f, Delta] = deal(max(Delta_f, Delta_z));
      end
      
      % Print the details
      deft_funnel_printout(nit, evaluations, iterate, setting,              ...
          indicators, vmax, norm_d_x, norm_d_s, Delta_f, Delta_z, rho,      ...
          it_type, '', sample_set);
      
      % Start a new iteration
      continue

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%      NORMAL STEP      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [nstep, nstep_x, nstep_s, exit_nstep] =                                 ...
        deft_funnel_normal_step(iterate, setting, indicators, models, Delta_z);

    if (exit_nstep == 2) % An infeasible stationary point has been found
        return
    end
    norm_n  = norm(nstep, Inf);
    
    % Compute the model of the objective function of the tangent problem
    M       = [models.derivatives.HLag zeros_nm; zeros_mn zeros_mm];
    gfxmod  = [models.derivatives.gfx; zeros_m1 ];
    m_xpn   = (gfxmod.' * nstep) + 0.5 * nstep.' * (M * nstep);
    g_n     = gfxmod + M * nstep;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%      TANGENT STEP      %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [norm_tstep, tstep, tstep_x, tstep_s, d, d_x, d_s, iterate,             ...
        indicators] = deft_funnel_tangent_step(iterate, setting,            ...
        indicators, models, Delta, nstep, nstep_x, nstep_s, M, g_n);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%      BUILD TRIAL POINT      %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    norm_d   = norm(d, Inf);
    norm_d_x = norm(d_x, Inf);
    norm_d_s = norm(d_s, Inf);

    % Compute the trial point
    iterate_plus = iterate;
    iterate_plus.x = iterate.x + d_x;
    iterate_plus.s = iterate.s + d_s;

    % Include the trial point in X and evaluate the functions f and c  
    % (xstatus(totalpoints) is set to 0 but is updated later on)
    if (any(d))
    
        % Update X and evaluate function at iterate_plus
        [sample_set, iterate_plus, evaluations, xstatus, sstatus] =         ...
            deft_funnel_eval_functions(f, c, h, iterate_plus, sample_set,   ...
            setting, evaluations, xstatus, const.unused, sstatus);

        % Compute the violation at the trial point
        iterate_plus.z_s = iterate_plus.zeval - iterate_plus.s;        
        vplus = 0.5*(iterate_plus.z_s.' * iterate_plus.z_s);
         
    end
  
    % Compute the model of f at the trial point
    m_xpd = gfxmod.' * d + 0.5 * d.' * (M * d);

    % Compute the improvement in the model of f at the trial point
    delta_f = -m_xpd;

    % Compute the improvement in the model of f at the trial point 
    % over that after the normal step
    delta_ft = m_xpn - m_xpd;
    
    if (delta_ft < 0)
        norm_tstep
        fprintf(' Warning: delta_ft = %.4e\n', delta_ft);
    end
   
    rho = (iterate.feval-iterate_plus.feval+setting.rho_eps)/(delta_f+setting.rho_eps);
    
    if (rho >= setting.eta1)
        succ = 1;
    else
        succ = 0;
    end
    
    % Compute the improvement in the model of c at the trial point
    Jsd = iterate.z_s + models.derivatives.J_s * d;
    Jsn = iterate.z_s + models.derivatives.J_s * nstep;
    
    delta_vsq  = indicators.v - 0.5 * (Jsd.' * Jsd);
    delta_vsqn = indicators.v - 0.5 * (Jsn.' * Jsn);
    
    pos = 0; % indicates if iterate_plus has been included in Y

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%      MU - ITERATION      %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (norm_d == 0.0)
        
        it_type = 'mu-iteration';
        rho     = 0;
        mu_iter = mu_iter + 1;
        
        if (setting.verbose >= 2)
            disp(' Mu-iteration:')
            disp(' ')
        end

        [sample_set, evaluations, Delta_f, Delta_z, xstatus, sstatus,       ...
            it_type, poised_model] = deft_funnel_improve_poisedness(f, c,   ...
            h, sample_set, iterate, setting, evaluations, indicators,       ...
            const, xstatus, sstatus, exit_nstep, Delta_f, Delta_z,          ...
            it_type, poised_model);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%      F - ITERATION      %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    elseif (norm_tstep > 0.0 && delta_f >= setting.kappa_delta * delta_ft && ...
            vplus <= vmax)
        
        mu_iter = 0;
        
        % =====================================================================
        % ===============      F - Successful iteration      ==================
        % =====================================================================
        
        if (rho >= setting.eta1)
            
            it_type = 'f-succ' ;
            
            if (setting.verbose >=2)
                disp(' Successful f-iteration')
                disp(' ')
            end

            [sample_set, iterate, setting, pos, xstatus, Delta_f,           ...
                Delta_z, Delta_f_counter, vmax, poised_model] =             ...
                deft_funnel_succ_iteration(1, sample_set, iterate,          ...
                iterate_plus, setting, indicators, model_size, Delta_f,     ...
                Delta_z, Delta_f_counter, Deltamax, vplus, vmax,            ...
                xstatus, const, succ, norm_d, norm_d, rho, poised_model);

        else
            
            % =================================================================
            % ==============      Second-order correction      ================
            % =================================================================
            
            % Compute the 2nd-order correction as the minimizer of 
            % || z(x+s) + J(x,s) n || subject to || n || <= Delta_z and
            % other constraints
            % (for large problems this should be done approximately)
            
            try_soc = 1 ;
            
            if (try_soc == 1)
                
                if (setting.verbose >=2)
                    disp(' First trial point was unsuccessful.')
                    disp(' Computing 2nd-order correction.')
                    disp(' ')
                end
                
                [sample_set, iterate, setting, evaluations, xstatus,        ...
                    sstatus, Delta_f, Delta_z, Delta_f_counter, vmax, rho,  ...
                    pos, it_type, poised_model] =                           ...
                    deft_funnel_sec_order_correction(f, c, h, sample_set,   ...
                    iterate, iterate_plus, setting, models, indicators,     ...
                    evaluations, d, delta_f, vmax, pos, xstatus, sstatus,   ...
                    const, model_size, Delta, Delta_f, Delta_z, Deltamax,   ...
                    Delta_f_counter, rho, succ, it_type, poised_model);
            end
        end
        
        % =====================================================================
        % ==============      F - Unsuccessful iteration      =================
        % =====================================================================
        
        % Enter if iteration is unsuccessful or the point couldn't 
        % be included in Y yet
        if (rho < setting.eta1 || pos == 0)
            
            it_type = 'f-unsucc' ;
            
            if (setting.verbose >= 2 && rho < setting.eta1)
                disp(' Unsuccessful f-iteration')
                disp([' rho = ', num2str(rho)])
                disp([' iterate.feval = ', num2str(iterate.feval)])
                disp([' iterate_plus.feval = ', num2str(iterate_plus.feval)])
                disp([' delta_f = ', num2str(delta_f)])
                disp([' delta_ft = ', num2str(delta_ft)])
                disp(' ')
            end
            
            [sample_set, iterate, setting, vmax, xstatus, it_type,          ...
                Delta_f, Delta_f_counter, Delta_z, poised_model,            ...
                changedPoint, newPoint, distNewPoint] =                     ...
                deft_funnel_unsucc_iteration(pos, 1, sample_set,            ...
                iterate, iterate_plus, setting, indicators, model_size,     ...
                const, xstatus, Delta, Delta_f, Delta_z, Delta_f_counter,   ...
                Delta_max_counter, Deltamax, delta_vsq, delta_vsqn,         ...
                vplus, vmax, rho, succ, norm_d, norm_n, it_type,            ...
                poised_model);

        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%      Z - ITERATION      %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    else
        
        mu_iter = 0;
        
        if (setting.verbose > 1)
            if (norm_tstep == 0.0)
                fprintf('  **** z-iteration because norm(t) = 0.0 \n');
            end
            if (delta_f < setting.kappa_delta * delta_ft)
                fprintf('%s %8.2e %s %8.2e\n',                              ...
                        '  **** z-iteration because delta_f = ',            ...
                        delta_f, ' < frac * delta_ft = ', delta_ft);
            end
            if (vplus > vmax)
                fprintf('%s %8.2e %s %8.2e\n',                              ...
                        '  **** z-iteration because vplus =', vplus,        ...
                        ' > vmax = ', vmax);
            end
        end
        
        rho = (indicators.v-vplus+setting.rho_eps)/(delta_vsq + setting.rho_eps);
        
        if (indicators.v-vplus < setting.rho_eps && delta_vsq < setting.rho_eps)
            rho = -1;
        end
        
        % =====================================================================
        % ===============      Z - Successful iteration      ==================
        % =====================================================================
        
        if (rho >= setting.eta1 && delta_vsq >= 0.4 * delta_vsqn && norm_n ~= 0)
            
            it_type = 'z-succ';
            succ    = 1;
            
            if (setting.verbose >=2)
                disp(' Successful z-iteration')
                disp(' ')
            end
            
            [sample_set, iterate, setting, pos, xstatus, Delta_z,           ...
                ~, Delta_z_counter, vmax, poised_model] =                   ...
                deft_funnel_succ_iteration(2, sample_set, iterate,          ...
                iterate_plus, setting, indicators, model_size, Delta_z,     ...
                Delta_z, Delta_z_counter, Deltamax, vplus, vmax,            ...
                xstatus, const, succ, norm_n, norm_d, rho, poised_model);
        end
        
        % =====================================================================
        % ==============      Z - Unsuccessful iteration      =================
        % =====================================================================
        
        % Enter if the iteration is unsuccessful or the point couldn't be 
        % included in Y yet
        if ((rho < setting.eta1 || delta_vsq < 0.2*delta_vsqn || pos == 0) && ...
               ~strcmp(it_type, 'mu-iteration'))
            
            it_type = 'z-unsucc';
            
            if (setting.verbose >= 2)
                disp(' Unsuccessful z-iteration')
                disp(' ')
            end
            
            if(rho < setting.eta1 || delta_vsq < 0.2 * delta_vsqn)
                succ = 0;
            else
                succ = 1;
            end
            
            [sample_set, iterate, setting, vmax, xstatus, it_type,          ...
                Delta_z, Delta_z_counter, Delta_z, poised_model,            ...
                changedPoint, newPoint, distNewPoint] =                     ...
                deft_funnel_unsucc_iteration(pos, 2, sample_set,            ...
                iterate, iterate_plus, setting, indicators, model_size,     ...
                const, xstatus, Delta, Delta_z, Delta_z, Delta_z_counter,   ...
                Delta_max_counter, Deltamax, delta_vsq, delta_vsqn,         ...
                vplus, vmax, rho, succ, norm_n, norm_n, it_type,            ...
                poised_model);
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%      UPDATE MODELS      %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    [models, iterate, indicators] =                                         ...
        deft_funnel_update_iteration(sample_set, iterate, dev_f, dev_h,     ...
        setting);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%      Iteration printout      %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sample_set = deft_funnel_poisedness_Y(sample_set, iterate, setting);
    
    deft_funnel_printout(nit, evaluations, iterate, setting,                ...
        indicators, vmax, norm_d_x, norm_d_s, Delta_f, Delta_z, rho,        ...
        it_type, '', sample_set);
    
    % Test for non-suitable model
    if (strcmp(setting.type_f, 'BB'))
        if (~isempty(find(isnan(models.f), 1))                           || ...
            ~isempty(find(~isreal(models.f), 1))                         || ...
            ~isempty(find(isinf(models.f), 1)))
            
            msg = 'Error: model_f contains NaN or Inf or nonreal components!!';
            if (setting.verbose >= 1)
                disp(msg);
            end
            return
        end
    end
    if (setting.cons_c)
        if (isa(c, 'function_handle'))
            if (~isempty(find(isnan(models.c), 1))                           || ...
                ~isempty(find(~isreal(models.c), 1))                         || ...
                ~isempty(find(isinf(models.c), 1)))

                msg = 'Error: model_c contains NaN or Inf or nonreal components!!';
                if (setting.verbose >= 1)
                    disp(msg);
                end
                return
            end
        end
    end
    
    % Test for the maximum number of function evaluations.
    if (evaluations.nfeval >= setting.maxeval)
       if (setting.verbose >= 1)
        	disp([' !!!!!!!!!!!!!!!!!!!!!!!!!!!!! MAX. EVALUATIONS',        ...
                  ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
        end
        return
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %%%%%%%%%%%%%%%%%      END OF MAIN ITERATION      %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (setting.verbose >= 1)
   disp([' !!!!!!!!!!!!!!!!!!!!!!!!!!!!! MAX. ITERATIONS',                  ...
         ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %%%%%%%%%%%%%%%      END OF DEFT_FUNNEL_MAIN      %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
