function [sample_set, setting, xstatus, pos] =                              ...
    deft_funnel_update_succ_iter_Y(sample_set, iterate_plus, setting,       ...
    model_size, xstatus, const, succ)
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Desc: Tries to include the successful trial point into the interpolation set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Augment interpolation set if not fully quadratic yet
if (setting.cur_degree < model_size.pquad || (setting.whichmodel == 3 &&    ...
    setting.cur_degree < 2*model_size.pquad))
    
    [sample_set, pY] = deft_funnel_augment_Y(sample_set, iterate_plus.x, setting);
    setting.cur_degree = pY;
    pos = pY;
    
else

    % Include xplus in the interpolation set, by replacing
    % another point if the model is already fully quadratic.
    [sample_set, pos] = deft_funnel_include_in_Y(sample_set,                ...
        iterate_plus.x, setting, [1:setting.cur_degree], setting.Lambda_XN, ...
        setting.criterion_S, succ);
    
    if (pos > 0)
        xstatus(sample_set.ind_Y(pos))  = const.unused;
    end
    
end

if (pos > 0)
    xstatus(sample_set.nbPoints) = const.inY;
    sample_set.ind_Y(pos) = sample_set.nbPoints;
    if (strcmp(setting.type_f, 'BB'))
        sample_set.fY(pos) = iterate_plus.feval;
    end
    if (setting.cons_c)
        sample_set.cY(:, pos) = iterate_plus.ceval;
    end
    if (setting.verbose >= 2)
        disp([' Replacement/inclusion of interpolation point at position ', ...
            int2str(pos), ' in Y successful'])
        disp(' ')
    end
end
    
end % end of deft_funnel_update_succ_iter_Y
