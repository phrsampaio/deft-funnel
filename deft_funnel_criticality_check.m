function critical = deft_funnel_criticality_check( iterate, indicators,     ...
    setting, it_type, Delta, norm_d, level )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Checks if convergence is attained in the full space or subspace.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

factor_opt = 1.0e+0;
if ( strcmp( level, 'subspace' ) )
    if ( ( indicators.norm_c_s <= factor_opt*setting.epsilon && indicators.norm_glag <= factor_opt*setting.epsilon ) || ...
          Delta <= setting.epsilon * setting.stallfact_subsp * norm( iterate.x ) )
        critical = 2;
    else
        critical = 0;
    end

else
    if ( ( indicators.norm_c_s <= setting.epsilon_i && indicators.norm_glag <= setting.epsilon_i ) || ...
         ( ~strcmp( it_type, '' )  && ~strcmp( it_type, 'mu-iteration' ) && ~strcmp( it_type, 'mu-iteration: improving poisedness' ) && ...
            ( Delta <= setting.epsilon * setting.stallfact_fullsp * norm( iterate.x ) || ...
              ( norm_d <= setting.epsilon * setting.stallfact_fullsp * norm( iterate.x ) && ~strcmp( it_type, 'criticality test' ) ) ) ) )

        critical = 1;
    else
        critical = 0;
    end 
end
           
end % end of deft_funnel_criticality_check