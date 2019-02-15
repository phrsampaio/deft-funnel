function errg = deft_funnel_compute_error( sampleSet, iterate, setting )
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Computes the bound for the error betweeen the gradient of the true 
% functions and the gradient of the surrogate models.
%
% Input: 
%   - sampleSet : struct of the sample set
%   - iterate   : struct of the current iterate
%   - setting   : struct of parameters
%
% Output:
%   - errg      : estimated error between the gradients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
sampleSet = deft_funnel_poisedness_Y( sampleSet, iterate, setting );
errg = sampleSet.lambda * sampleSet.Y_radius;
                
end % end of deft_funnel_compute_error