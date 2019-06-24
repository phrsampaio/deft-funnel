function sampleSet = deft_funnel_replace_in_Y( sampleSet, ynew, j, setting )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Updates the interpolation set for a transformation of Y in Yplus where 
% the vector ynew replaces Y(:,j).  Also updates the factorization of Z(Y)
% accordingly.
%
% Input:
%   - sampleSet : struct of the sample set
%   - ynew      : the vector hich is to replace Y(:,j) in the interpolation set
%   - j         : the index of the interpolation point to be replaced
%   - setting   : struct of parameters
%
% Output:
%   - sampleSet : updated struct of the sample set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Replace new point in Y
sampleSet.Y( :, j ) = ynew;

% Compute new factorization
sampleSet = deft_funnel_build_QR_of_Y( sampleSet, setting );

end % end of deft_funnel_replace_in_Y

