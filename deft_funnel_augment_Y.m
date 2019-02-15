function [ sampleSet, pY ] = deft_funnel_augment_Y( sampleSet, Ynew, setting )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Augments the interpolation set by adding new point(s). This assumes 
% that the model is not yet fully quadratic. If this is the case,
% the current interpolation set (and the associated polynomial degree) are
% increased by the number of columns in Ynew.
%
% Input:
%   - sampleSet : struct of the sample set
%   - Ynew      : new point(s)
%   - setting   : struct of parameters
%
% Output:
%   - sampleSet : updated struct of the sample set
%   - pY        : update number of sample points in Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ n, pY ] = size( sampleSet.Y );

if ( ( pY >= ( ( n + 1 ) * ( n + 2 ) ) / 2 ) && ( setting.whichmodel ~= 3 ) )

   disp( ' === augment_Y: warning!!! The interpolation is already fully quadratic!')
   disp( '     Ignoring augmentation...')

else

   sampleSet.Y = [ sampleSet.Y Ynew ];
   pY = pY + size( Ynew, 2 );
   
   sampleSet = deft_funnel_build_QR_of_Y( sampleSet, setting );

end

end % end of deft_funnel_augment_Y