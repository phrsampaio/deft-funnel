function [ sampleSet, sample, evaluations, xstatus, sstatus ] =             ...
           deft_funnel_augmX_evalfc( f, c, sample, sampleSet, setting,      ...
           evaluations, xstatus, xstatus_val, sstatus )
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Adds new 'sample' into the set of all points X. The objective 
% and constraint functions are evaluated at the new point and the 
% function values are added to fX and cX.
%
% Input:
%   - f            : objective function handle
%   - c            : constraints' function handle
%   - sample       : struct of the new point
%   - sampleSet    : struct of the sample set
%   - setting      : struct of parameters
%   - evaluations  : struct containing info about nb of evaluations
%   - xstatus      : indicates whether the points are in Y or not (1 or 0)
%   - xstatus_val  : value (1 or 0) associated to the new point in 'xstatus'
%   - sstatus      : indicates whether the points are in the subspace 
%                    or not (1 or 0)
%
% Output:
%   - sampleSet    : updated struct of the sample set
%   - sample       : updated with the function values evaluated
%   - evaluations  : updated based on the new point
%   - xstatus      : updated based on the new point
%   - sstatus      : updated based on the new point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

full_n  = sample.fulldim;
I       = eye(full_n);
xfix    = sample.xfix;
nfix    = sample.nfix;
indfix  = sample.indfix;
indfree = sample.indfree;

% Update the number of evaluated points and set the index of the new point
sampleSet.nbPoints = sampleSet.nbPoints + 1;
index = sampleSet.nbPoints;

% Set status of the new point
xstatus(index) = xstatus_val;
sstatus(index) = 1; % the point is contained in current subspace

% Augment X with full-dimensional y, scale if user-defined and 
% evaluate f and c at y
if ( nfix > 0 )

    samplefull  = I(:,indfix) * xfix(indfix) + I(:,indfree) * sample.x;
    sampleSet.X(:,index) = samplefull;

    if ( setting.scaleX )
      samplefull = samplefull ./ setting.scalefacX;
    end

    try
        fvalue = f(samplefull);
    catch
        disp(' ')
        disp( ' Error: evaluation of objective function FAILED at the point');
        samplefull
        disp( ' Possible causes:' )
        disp([ ' (1) Hidden constraint: your black-box function might not be able to be evaluated', ...
               ' at that point.' ] );
        disp([ ' (2) Function handle: check if the function handle', ...
               ' that you have passed as input is correct.'] );
        disp([ ' (3) Black-box function: check if your black-box', ...
               ' function does not contain any internal errors.'] );
        disp( ' Setting fvalue = 1.0e+10. This might affect the final solution.' );
        disp(' ')
        fvalue = 1.0e+10;
    end

    try
        cvalue = c(samplefull);
    catch
        disp(' ')
        disp(' Error: evaluation of constraint function FAILED at the point' );
        samplefull
        disp( ' Possible causes:' )
        disp([ ' (1) Hidden constraint: your black-box function might not be able to be evaluated', ...
               ' at that point.' ] );
        disp([ ' (2) Function handle: check if the function handle', ...
               ' that you have passed as input is correct.'] );
        disp([ ' (3) Black-box function: check if your black-box', ...
               ' function does not contain any internal errors.'] );
        disp( ' Setting cvalue = 1.0e+10*ones(m,1). This might affect the final solution.' );
        disp(' ')
        cvalue = 1.0e+10 * ones( sample.sdim, 1 );
    end

else

    sampleSet.X(:,index) = sample.x;
    
    if ( setting.scaleX )
        sample.x = sample.x ./ setting.scalefacX;
    end

    try
        fvalue = f(sample.x);
    catch
        disp(' ')
        disp(' Error: evaluation of objective function FAILED at the point');
        sample.x
        disp( ' Possible causes:' )
        disp([ ' (1) Hidden constraint: your black-box function might not be able to be evaluated', ...
               ' at that point.' ] );
        disp([ ' (2) Function handle: check if the function handle', ...
               ' that you have passed as input is correct.'] );
        disp([ ' (3) Black-box function: check if your black-box', ...
               ' function does not contain any internal errors.'] );
        disp(' Setting fvalue = 1.0e+10. This might affect the final solution.');
        disp(' ')
        fvalue = 1.0e+10;
    end

    try
        cvalue = c(sample.x);
    catch
        disp(' ')
        disp(' Error: evaluation of constraint function FAILED at the point');
        sample.x
        disp( ' Possible causes:' )
        disp([ ' (1) Hidden constraint: your black-box function might not be able to be evaluated', ...
               ' at that point.' ] );
        disp([ ' (2) Function handle: check if the function handle', ...
               ' that you have passed as input is correct.'] );
        disp([ ' (3) Black-box function: check if your black-box', ...
               ' function does not contain any internal errors.'] );
        disp(' Setting cvalue = 1.0e+10*ones(m,1). This might affect the final solution.');
        disp(' ')
        cvalue = 1.0e+10 * ones( sample.sdim, 1 );
    end

end

% Augment fX with new function value
if ( max( size( fvalue ) ) > 1 )
   sampleSet.fX(index) = min( setting.fmax, norm( fvalue )^2);
else
   sampleSet.fX(index) = min( setting.fmax, real( fvalue ) );
end
sample.feval = sampleSet.fX(index);

% Augment cX with new function value
length_c = max( size( cvalue ) );
if ( length_c > 1 )
   for i=1:length_c
     sampleSet.cX(i,index) = min( setting.fmax, real( cvalue(i) ) );
     sample.ceval(i) = sampleSet.cX(i,index);
   end
else
   sampleSet.cX(index) = min( setting.fmax, real( cvalue ) );
   sample.ceval = sampleSet.cX(index);
end

% Update the number of evaluations and points evaluated
evaluations.nfeval = evaluations.nfeval + 1;
evaluations.nceval = evaluations.nceval + 1;

end % end of deft_funnel_augmX_evalfc
