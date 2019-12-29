function models = deft_funnel_build_models(sample_set, setting)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Builds the interpolating models for the objective and the constraint
% functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = size(sample_set.cY, 1);

% Check if the objective function is of black-box type
if (strcmp(setting.type_f, 'BB'))
    model_f  = deft_funnel_computeP(sample_set, sample_set.fY, setting);
    models.f = model_f;
end

% Check if there are black-box constraints
if (setting.cons_c)
    for i=1:m
        model_c(i,:) = deft_funnel_computeP(sample_set, sample_set.cY(i,:), ...
            setting);
        models.c = model_c;
    end
end
          
end % end of deft_funnel_build_models
