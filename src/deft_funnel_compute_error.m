function errg = deft_funnel_compute_error(sample_set, iterate, setting)
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Computes the bound for the error betweeen the gradient of the true 
% functions and the gradient of the surrogate models.
%
% Input: 
%   - sample_set : struct of the sample set
%   - iterate    : struct of the current iterate
%   - setting    : struct of parameters
%
% Output:
%   - errg      : estimated error between the gradients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample_set = deft_funnel_poisedness_Y(sample_set, iterate, setting);
errg = sample_set.lambda * sample_set.Y_radius;
                
end % end of deft_funnel_compute_error