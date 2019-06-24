function sampleSet = deft_funnel_swap_in_Y( i, j, sampleSet )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Swaps the position of interpolation points i and j in Y, and updates 
% the factorization of Z(Y) accordingly.
%
% Input:
%   - i         : the position in which Y(:,j) should be inserted
%   - j         : the position in which Y(:,i) should be inserted
%   - sampleSet : sample set structure
%
% Output:
%   - sampleSet : sample set structure updated
%
% Called by: deft_funnel_succ_iteration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensure that ii is smaller than jj.
if ( i > j )
   ii = j;
   jj = i;
elseif ( i < j )
   ii = i;
   jj = j;
else
   return
end

% Permute the columns of Y, the indices in ind_Y and the values in fY and cY
y                   = sampleSet.Y(:,ii);
sampleSet.Y(:,ii)   = sampleSet.Y(:,jj);
sampleSet.Y(:,jj)   = y;

ind                 = sampleSet.ind_Y(ii);
sampleSet.ind_Y(ii) = sampleSet.ind_Y(jj);
sampleSet.ind_Y(jj) = ind;

f                   = sampleSet.fY(ii);
sampleSet.fY(ii)    = sampleSet.fY(jj);
sampleSet.fY(jj)    = f;

c                   = sampleSet.cY(:,ii);
sampleSet.cY(:,ii)  = sampleSet.cY(:,jj);
sampleSet.cY(:,jj)  = c;

end % end of deft_funnel_swap_in_Y