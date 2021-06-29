function sample_set = deft_funnel_swap_in_Y(i, j, sample_set, setting)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Swaps the position of interpolation points i and j in Y, and updates 
% the factorization of Z(Y) accordingly.
%
% Input:
%   - i          : the position in which Y(:,j) should be inserted
%   - j          : the position in which Y(:,i) should be inserted
%   - sample_set : struct of the sample set
%   - setting    : struct of parameters
%
% Output:
%   - sample_set : struct of the sample set updated
%
% Called by: deft_funnel_succ_iteration, deft_funnel_unsucc_iteration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensure that ii is smaller than jj.
if (i > j)
   ii = j;
   jj = i;
elseif (i < j)
   ii = i;
   jj = j;
else
   return
end

% Permute the columns of Y, the indices in ind_Y and the values in fY and cY
y                    = sample_set.Y(:,ii);
sample_set.Y(:,ii)   = sample_set.Y(:,jj);
sample_set.Y(:,jj)   = y;

ind                  = sample_set.ind_Y(ii);
sample_set.ind_Y(ii) = sample_set.ind_Y(jj);
sample_set.ind_Y(jj) = ind;

if (strcmp(setting.type_f, 'BB'))
    f                 = sample_set.fY(ii);
    sample_set.fY(ii) = sample_set.fY(jj);
    sample_set.fY(jj) = f;
end

if (setting.cons_c)
    c                   = sample_set.cY(:,ii);
    sample_set.cY(:,ii) = sample_set.cY(:,jj);
    sample_set.cY(:,jj) = c;
end

end % end of deft_funnel_swap_in_Y
