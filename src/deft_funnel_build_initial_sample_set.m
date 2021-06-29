function [sample_set, iterate, evaluations, xstatus, sstatus, poised_model, ...
    msg] = deft_funnel_build_initial_sample_set(f, c, h, sample_set,        ...
    iterate, setting, evaluations, xstatus, sstatus, const, model_size,     ...
    Delta)
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Builds the initial sample set using either of the followinbg 
% two stratagies:
% (1) random interpolation points
% (2) simplex approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (setting.verbose >= 2)
    disp(' *********** deft_funnel_build_initial_sample_set ***********')
    disp(' *** Building initial sample set')
    disp([ ' *** Degree of the initial model = ', int2str(setting.cur_degree)])
end

% Compute an initial poised interpolation set around the starting point
sample_set.Y(:, 1) = iterate.x;
if (strcmp(setting.initial_Y, 'random'))

    sample_set.Y(:, 2:setting.cur_degree) =                                 ...
        -ones(iterate.xdim, setting.cur_degree - 1) +                       ...
        2*rand(iterate.xdim, setting.cur_degree - 1);

    for  j = 2:setting.cur_degree
        sample_set.Y(:, j) = sample_set.Y(:, 1) +                           ...
            sample_set.Y(:, j) * (Delta / norm(sample_set.Y(:, j)));

        % Make sure that the interpolation points are within the bounds
        sample_set.Y(:, j) = max(min(sample_set.Y(:, j), setting.ux), setting.lx);
    end

    % Build the initial factorization
    sample_set = deft_funnel_build_QR_of_Y(sample_set, setting);

    % Make the interpolation set poised
    sample_set = deft_funnel_repair_Y(sample_set, iterate, setting, Delta);

elseif strcmp(setting.initial_Y, 'simplex')

    % Compute the initial interpolation set (simplex plus midpoints)
    I = eye(iterate.xdim);

    for j = 1:iterate.xdim

        % Linear
        step1 = -Delta;
        sample_set.Y(:, j + 1) = iterate.x + step1 * I(:, j);

        % Diagonal
        if (setting.cur_degree >= model_size.pdiag)
            step2 =  Delta;
            sample_set.Y(:, j + 1 + iterate.xdim) = iterate.x + step2 * I(:, j);
        end
    end

    % Quadratic
    if (setting.cur_degree == model_size.pquad)
        k = 2 * iterate.xdim + 2;
        for j = 1:iterate.xdim
            for jj = j+1:iterate.xdim
                sample_set.Y(:, k) = 0.5 * (sample_set.Y(:, j + 1) +        ...
                    sample_set.Y(:, jj + 1));
                k = k + 1;
            end
        end
    end

    % Make sure that the interpolation points are within the bounds
    for  j = 2:setting.cur_degree
        sample_set.Y(:, j) = max(min(sample_set.Y(:, j), setting.ux), setting.lx);
    end
   
    % Build the initial factorization
    sample_set = deft_funnel_build_QR_of_Y(sample_set, setting);
end 

% Compute the models gradient error associated to the sample set by
% computing the poisedness of the sample set
sample_set = deft_funnel_poisedness_Y(sample_set, iterate, setting);

% Compute the associated function values, possibly reducing Delta to
% ensure that the objective function remains finite at all interpolation
% points.

% Set a sample structure
sample = deft_funnel_create_sample([], [], iterate);

for i = 1:setting.cur_degree

   % Store the new points in X and evaluate f and c at them
   sample.x = sample_set.Y(:, i);
   [sample_set, sample, evaluations, xstatus, sstatus] =                  	...
       deft_funnel_eval_functions(f, c, h, sample, sample_set, setting,     ...
       evaluations, xstatus, const.inY, sstatus);

   % Retrieve function values of the iterate
   if (i == 1)
       iterate.feval = sample.feval;
       iterate.ceval = sample.ceval;
       iterate.zeval = sample.zeval;
   end
   
   if (strcmp(setting.type_f, 'BB' ))
       sample_set.fY(i)    = sample_set.fX(i);
   end
   
   if (setting.cons_c)
       sample_set.cY(:, i) = sample_set.cX(:, i);
   end

   if (evaluations.nfeval > setting.maxeval || evaluations.nceval > setting.maxeval)
       msg = [' Error: Maximum number of ', int2str(setting.maxeval),       ...
              ' function evaluations reached.'];
       % Including fixed variables at return
       if (iterate.nfix > 0)
           I = eye(iterate.xdim + iterate.nfix);
           iterate.x = I(:, iterate.indfix) * setting.lx(iterate.indfix) +  ...
               I(:, iterate.indfree) * iterate.x;
       end
       disp(msg)
       poised_model = 0;
       return
   end

end

xstatus = xstatus';
poised_model = 1;

sample_set.ind_Y = 1:setting.cur_degree;
sample_set.i_xbest = sample_set.ind_Y(1);

msg = 'Initial sample set successfully built';

if (setting.verbose >= 2)
    disp(' *** Initial sample set successfully built')
    disp(' ***** return from deft_funnel_build_initial_sample_set *****')
end

end % end of deft_funnel_build_initial_sample_set
