function [sample_set, sample, evaluations, xstatus, sstatus] =              ...
    deft_funnel_eval_functions(f, c, h, sample, sample_set, setting,        ...
    evaluations, xstatus, xstatus_val, sstatus, varargin)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Evaluates the objective and constraint functions at 'sample', 
% adds the function values to fX and cX, and adds the new 'sample' into 
% the set of all points X.
%
% Input:
%   - f            : objective function handle
%   - c            : BB constraints function handle or string "combined" that
%                    indicates that the black box returns both objective and 
%                    contraint functions simultaneously
%   - h            : function handle to the white-box constraints
%   - sample       : struct of the new point
%   - sample_set   : struct of the sample set
%   - setting      : struct of parameters
%   - evaluations  : struct containing info about nb of evaluations
%   - xstatus      : indicates whether the points are in Y or not (1 or 0)
%   - xstatus_val  : value (1 or 0) associated to the new point in 'xstatus'
%   - sstatus      : indicates whether the points are in the subspace 
%                    or not (1 or 0)
%
% Output:
%   - sample_set   : updated with the new sample point
%   - sample       : updated with the evaluations of the functions f and c
%   - evaluations  : updated with new function evaluations
%   - xstatus      : updated with status of the new points
%   - sstatus      : updated with status of the new points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eval_only = false;
if (nargin > 10)
    if (strcmp(varargin{1}, 'eval_only'))
        eval_only = true;
    end
end

full_n  = sample.fulldim;
I       = eye(full_n);
xfix    = sample.xfix;
nfix    = sample.nfix;
indfix  = sample.indfix;
indfree = sample.indfree;

if (~eval_only)
    % Update the number of evaluated points and set the index of the new point
    sample_set.nbPoints = sample_set.nbPoints + 1;
    index = sample_set.nbPoints;

    % Set status of the new point
    xstatus(index) = xstatus_val;
    sstatus(index) = 1; % the point is contained in the current subspace
end

% Error message in case of black-box evaluation failure
error_msg = [' Possible causes:\n',                                         ...
             ' (1) Hidden constraint: your function might not',             ...
             ' be able to be evaluated at that point.\n',                   ...
             ' (2) Function handle: check if the function handle',          ...
             ' that you have passed as input is correct.\n',                ...
             ' (3) Black-box function: check if your black-box',            ...
             ' function does not contain any internal errors.\n',           ...
             ' Setting output(s) to 1.0e+10. This might affect the',        ...
             ' final solution.\n\n'];
         
big_value = 1.0e+10; % Used as function value when the evaluation fails

% Augment X with full-dimensional iterate, scale it if user-defined and 
% evaluate f and c at it
if (nfix > 0)
    my_sample = I(:,indfix) * xfix(indfix) + I(:,indfree) * sample.x;
else
    my_sample = sample.x;
end

if (~eval_only)
    sample_set.X(:,index) = my_sample;
end

if (setting.scaleX)
  my_sample = my_sample ./ setting.scalefacX;
end

if (strcmp(c, 'combined')) % Single black-box call for both f and c

    try
        output = f(my_sample);
        fvalue = output(1);
        cvalue = output(2:setting.nb_cons_c+1);
    catch
        disp(' ')
        disp( ' Error: evaluation of the black box FAILED at the point');
        mysample
        fprintf(error_msg);
        fvalue = big_value;
        cvalue = big_value * ones(setting.nb_cons_c, 1);
    end

else

    try
        fvalue = f(my_sample);
    catch
        disp(' ')
        disp(' Error: evaluation of objective function FAILED at the point' );
        my_sample
        fprintf(error_msg);
        fvalue = big_value;
    end

    % Check if there are black-box constraints
    if (setting.cons_c)
        try
            cvalue = c(my_sample);
        catch
            disp(' ')
            disp(' Error: evaluation of constraint function c FAILED at the point' );
            my_sample
            fprintf(error_msg);
            cvalue = big_value * ones(setting.nb_cons_c, 1);
        end
    else
        cvalue = [];
    end

end

% Check if there are white-box constraints
if (setting.cons_h)
    try
        hvalue = h(my_sample);
    catch
        disp(' ')
        disp(' Error: evaluation of constraint function h FAILED at the point' );
        my_sample
        fprintf(error_msg);
        hvalue = big_value * ones(setting.nb_cons_h, 1);
    end
else
    hvalue = [];
end

% Store new f value and (possibly) add it to fX
if (max(size(fvalue)) > 1)
    fvalue = min(setting.fmax, norm(fvalue)^2);
else
    fvalue = min(setting.fmax, real(fvalue));
end
sample.feval = fvalue;

if (~eval_only)
    if (strcmp(setting.type_f, 'BB'))
        sample_set.fX(index) = fvalue;
    end
end

% Augment cX with new c value
if (setting.cons_c)
    length_c = max(size(cvalue));
else
    length_c = 0;
end
if (length_c > 1)
    for i=1:length_c
        if (~eval_only)
            sample_set.cX(i,index) = min(setting.fmax, real(cvalue(i)));
        end
        sample.ceval(i) = min(setting.fmax, real(cvalue(i)));
    end
elseif (length_c == 1)
    if (~eval_only)
        sample_set.cX(index) = min(setting.fmax, real(cvalue));
    end
    sample.ceval = min(setting.fmax, real(cvalue));
end

% Store new h value
if (setting.cons_h)
    length_h = max(size(hvalue));
else
    length_h = 0;
end
if (length_h > 1)
    for i=1:length_h
        sample.heval(i) = min(setting.fmax, real(hvalue(i)));
    end
elseif (length_h == 1)
    sample.heval = min(setting.fmax, real(hvalue));
end

if (length_c > 0 && length_h > 0)
    sample.zeval = [reshape(sample.ceval,[],1); reshape(sample.heval,[],1)];
elseif (length_c > 0)
    sample.zeval = reshape(sample.ceval,[],1);
else
    sample.zeval = reshape(sample.heval,[],1);
end

% Update the number of evaluations and points evaluated
evaluations.nfeval = evaluations.nfeval + 1;
if (setting.cons_c)
    evaluations.nceval = evaluations.nceval + 1;
end
if (setting.cons_h)
    evaluations.nheval = evaluations.nheval + 1;
end

end % end of deft_funnel_eval_functions
