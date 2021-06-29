function sample_set = deft_funnel_replace_in_Y(sample_set, ynew, j, setting)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Updates the interpolation set for a transformation of Y in Yplus where 
% the vector ynew replaces Y(:,j).  Also updates the factorization of Z(Y)
% accordingly.
%
% Input:
%   - sample_set : struct of the sample set
%   - ynew       : the vector hich is to replace Y(:,j) in the interpolation set
%   - j          : the index of the interpolation point to be replaced
%   - setting    : struct of parameters
%
% Output:
%   - sample_set : updated struct of the sample set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Replace new point in Y
sample_set.Y(:, j) = ynew;

% Compute new factorization
sample_set = deft_funnel_build_QR_of_Y(sample_set, setting);

end % end of deft_funnel_replace_in_Y

