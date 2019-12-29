function [sample_set, replaced, maximprove, badcond] =                      ...
    deft_funnel_repair_Y(sample_set, iterate, setting, Delta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Repairs the interpolation set Y by first replacing the interpolation 
% points at a distance larger than setting.factor_FPR*Delta from Y(:,1), and 
% then replaces interpolation points with the best possible point obtained by 
% taking the best optimal replacement, until this improvement no longer causes 
% a relative simplex volume increase by at least threshold. For conventions 
% on how polynomials are represented, see the documentation of evalZ.
%
% Input:
%   - sample_set : struct of the sample set
%   - iterate    : struct of the current point
%   - setting    : struct of parameters
%   - Delta      : trust-region radius
%
% Output:
%   - sample_set : updated struct of the sample set
%   - replaced   : indices of replaced points
%   - maximprove : maximum improvement value obtained
%   - badcond    : matrix conditioning (0 for well cond and 1 otherise)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_improve_loops = 20; % the max number of times every close point is 
                        % considered for improvement
badcond           = 0;  % interpolation matrix conditioning
                        % 0 = well conditioned; 1 = ill conditioned

if (setting.verbose > 1)
    verbose = 1;
else
    verbose = 0;
end

if (verbose)
    disp([ ' --- enter deft_funnel_repair_Y for Delta = ', num2str(Delta), ' ---' ]);
    disp(' ')
end

[~, p1] = size(sample_set.Y);
replaced  = [];

% Compute the distance of the current interpolation points to the base point.
d = zeros(1, p1);
for j=2:p1
    d(j) = norm(sample_set.Y(:,j)- sample_set.Y(:,1));
end
[dsorted, jsorted] = sort(d, 'descend');

if (verbose)
    disp(' Distance of current interpolation points to the iterate: ')
    d
    dsorted
    jsorted
end

% First replace the distant points.
if (verbose)
    disp(' Trying to replace far points first ')
    disp(' ')
end

jmax = 0;
for j = 1:p1
    if (dsorted(j) > setting.factor_FPR*(1 + setting.eps_L)*Delta)
        jmax = jsorted(j);
      
        if (setting.hardcons == 1)
            [y, improvement] = deft_funnel_find_new_yj_bc(sample_set,       ...
                iterate, jmax, setting, Delta);
        else
            [y, improvement] = deft_funnel_find_new_yj(sample_set, jmax,    ...
                setting, Delta);
        end
      
        if (verbose)
            disp([' ==> lambda = ' num2str(improvement) ' at j=' num2str(jmax)])
        end
      
        if (improvement > setting.Lambda_FP)
         
            sample_set = deft_funnel_replace_in_Y(sample_set, y, jmax, setting);
            replaced  = [replaced jmax];
            d(jmax) = norm(y - sample_set.Y(:,1)); 
        end
      
    else
        sample_set.Y_radius = dsorted(j);
        break
    end
   
end

if (verbose)
    oldLambda = sample_set.lambda;
    oldRadius = sample_set.Y_radius;
	 sample_set = deft_funnel_poisedness_Y(sample_set, iterate, setting);
	 disp(' Results of attempt to replace far points:')
    disp([' ==> Number of points replaced: ', num2str(length(replaced))])
    replaced
    disp([' ==> poisedness(Y) before replacement = ', num2str(oldLambda)])
    disp([' ==> poisedness(Y) after replacement = ', num2str(sample_set.lambda)])
    disp([' ==> sample_set.Y_radius before replacement  = ', num2str(oldRadius)])
    disp([' ==> sample_set.Y_radius after replacement  = ', num2str(sample_set.Y_radius)])
    disp(' ')
end

% Perform a loop over possible optimal improvements.
if (verbose)
    disp(' Trying to replace any point based on the maximum possible improvement')
    disp(' ')
end

for  k = 1:max_improve_loops

    maximprove = 0;

    % Loop on all the possible replacements to find the best one.
    for j = 2:p1

        if (setting.hardcons == 1)
            [y, improvement] = deft_funnel_find_new_yj_bc(sample_set,       ...
                iterate, j, setting, Delta);
        else
            [y, improvement] = deft_funnel_find_new_yj(sample_set, j,       ...
                setting, Delta);
        end

        % Remember the current polynomial value, index and replacement point
        % that is the best so far
        if (verbose)
            disp([' ==> j = ', int2str(j), ' improve = ', num2str(improvement)])
            y
        end
      
        if (improvement > maximprove)
            maximprove = improvement;
            jmax       = j;
            ymax       = y;
        end
    end

    % If no significant improvement was found, return after updating the 
    % interpolation radius.
    if (maximprove < setting.Lambda_CP || jmax == 0)
        sample_set.Y_radius = max(d);
        if ( verbose > 0 )
            disp(' No significant improvement was found')
            disp([' maximprove(small) = ', num2str(maximprove), ...
                ', jmax= ' int2str(jmax)])
        end
        return
    end

    % Perform the best replacement.
    if (verbose > 0)
        disp([' maximprove= ' num2str(maximprove) ', jmax= ' int2str(jmax)])
    end

    sample_set = deft_funnel_replace_in_Y(sample_set, ymax, jmax, setting);
    d(jmax) = norm(ymax - sample_set.Y(:,1));

    if (length(find(replaced == jmax)) == 0)
        replaced = [replaced jmax];
    end

end

% Recompute sample_set.Y_radius after the exchanges
sample_set.Y_radius = max(d); 

if (verbose)
    sample_set = deft_funnel_poisedness_Y(sample_set, iterate, setting);
    disp(' Final results of repairing Y:')
    disp([' ==> Number of points replaced: ', num2str(length(replaced))])
    replaced
    disp([' ==> poisedness(Y) = ', num2str(sample_set.lambda)])
    disp([' ==> sample_set.Y_radius  = ', num2str(sample_set.Y_radius)])
    disp(' --- exit deft_funnel_repair_Y ---')
end

end % end of deft_funnel_repair_Y
