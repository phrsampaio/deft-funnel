function models = deft_funnel_build_models( sampleSet, setting )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Builds the interpolating models for the objective and the constraint
% functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = size( sampleSet.cY, 1 );

model_f  = deft_funnel_computeP( sampleSet, sampleSet.fY, setting );
for i=1:m
    model_c(i,:) = deft_funnel_computeP( sampleSet, sampleSet.cY(i,:), ...
        setting);
end

models.f = model_f;
models.c = model_c;
          
end % end of deft_funnel_build_models
