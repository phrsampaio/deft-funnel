function [ Sreduced, mSreduced, asc_order ] = deft_funnel_multistart_best_samples(S, mS, kappa)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Builds a reduced sample set with (kappa * nbsamples) points from S.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbsamples = size(S, 2);
nbreduced = kappa * nbsamples;

[mSsorted, asc_order] = sort(mS);
mSreduced = mSsorted(1:nbreduced);

Sreduced = S(:,asc_order);
Sreduced = Sreduced(:,1:nbreduced);

end % end of deft_funnel_best_samples