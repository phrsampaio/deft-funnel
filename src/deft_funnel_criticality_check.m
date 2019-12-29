function critical = deft_funnel_criticality_check(iterate, indicators,      ...
    setting, it_type, Delta, norm_d, level)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Checks if convergence is attained in the full space or subspace.
%
% Output: 
%   - critical: 0 - iterate is not a critical point
%               1 - iterate is a critical point at the full space
%               2 - iterate is a critical point at a subspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
factor_opt = 1.0e+0;
if (strcmp(level, 'subspace'))
    if ((indicators.norm_z_s <= factor_opt*setting.epsilon && indicators.norm_glag <= factor_opt*setting.epsilon) || ...
        Delta <= setting.epsilon * setting.stallfact_subsp * norm(iterate.x))
        critical = 2;
    else
        critical = 0;
    end

else
    if ((indicators.norm_z_s <= setting.epsilon_i && indicators.norm_glag <= setting.epsilon_i) || ...
         (~strcmp(it_type, '')  && ~strcmp(it_type, 'mu-iteration' ) && ~strcmp(it_type, 'mu-iteration: improving poisedness') && ...
            (Delta <= setting.epsilon * setting.stallfact_fullsp * norm(iterate.x) || ...
              (norm_d <= setting.epsilon * setting.stallfact_fullsp * norm(iterate.x) && ~strcmp(it_type, 'criticality test')))))

        critical = 1;
    else
        critical = 0;
    end 
end
           
end % end of deft_funnel_criticality_check