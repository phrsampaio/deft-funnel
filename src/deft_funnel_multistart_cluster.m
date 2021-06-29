function output = deft_funnel_multistart_cluster(S, mS, index, radius)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Checks if there is a neighbor S(i) to S(index) with better merit 
% function value (1 if yes and 0 otherwise).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbsamples = size(S, 2);
for i = 1:nbsamples
    if (i ~= index)
        distSi = norm(S(:, i) - S(:, index));
        if (mS(i) < mS(index) && distSi < radius)
            output = 1;
            return
        end
    end
end
output = 0; % there is no better point in the neighborhood

end % end of deft_funnel_multistart_cluster