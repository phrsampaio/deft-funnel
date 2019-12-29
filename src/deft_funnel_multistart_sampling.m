function samples = deft_funnel_multistart_sampling(distribution, dimension, ...
    nbsamples, lx, ux)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Sample points from the uniform distribution.
%
% Input:
%   - distribution : distribution type to sample from
%   - dimension    : dimension of the sample points
%   - nbsamples    : number of points to be sampled
%   - lx           : lower bounds
%   - ux           : upper bounds
%
% Output:
%   - samples      : (nbsamples x dimension) matrix containing the sample
%                    points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (strcmp(distribution, 'uniform'))
    
    % Generate 'nbsamples' points between bounds 'lx' and 'ux'.
    % If the lower bound or the upper bound is not a real number, consider -10
    % and 10, respectively.
    for i=1:dimension
        if ~isinf(lx(i))
            a = lx(i);
        else
            a = -10;
        end
        if ~isinf(ux(i))
            b = ux(i);
        else
            b = max(0,a) + 10;
        end
        samples(i,:) = (b-a).*rand(1,nbsamples) + a;
    end
end

end % end of deft_funnel_multistart_sampling