function [sample_set, pY] = deft_funnel_augment_Y(sample_set, Ynew, setting)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Augments the interpolation set by adding new point(s). This assumes 
% that the model is not yet fully quadratic. If this is the case,
% the current interpolation set (and the associated polynomial degree) are
% increased by the number of columns in Ynew.
%
% Input:
%   - sample_set : struct of the sample set
%   - Ynew       : new point(s)
%   - setting    : struct of parameters
%
% Output:
%   - sample_set : updated struct of the sample set
%   - pY         : update number of sample points in Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n, pY] = size(sample_set.Y);

if ((pY >= ((n + 1) * (n + 2))/2) && (setting.whichmodel ~= 3))

   disp(' === augment_Y: warning!!! The interpolation is already fully quadratic!')
   disp('     Ignoring augmentation...')

else

   sample_set.Y = [sample_set.Y Ynew];
   pY = pY + size(Ynew, 2);
   
   sample_set = deft_funnel_build_QR_of_Y(sample_set, setting);

end

end % end of deft_funnel_augment_Y