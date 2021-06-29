function deft_funnel_printout(nit, evaluations, iterate, setting,           ...
    indicators, vmax, norm_d_x, norm_d_s, Delta_f, Delta_z, rho, it_type, msg, sample_set)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Prints the iteration summary. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isempty(sample_set))
    errg = sample_set.errg;
end
      
if (setting.verbose >= 1)
    if (nit > 0)
        if (setting.show_errg)
            fprintf('%5d  %5d  %+.14e  %.2e  %.2e  %.4e %.2e %.2e  %.2e  %.2e  %.2e  %+.2e  %s\n', ...
                nit, evaluations.nfeval, iterate.feval,                     ...
                indicators.norm_glag, indicators.norm_z_s, vmax, errg,      ...
                norm_d_x, norm_d_s, Delta_f, Delta_z, rho, it_type);
            if(~isempty(msg))
                disp(msg);
            end
        else
            fprintf('%5d  %5d  %+.14e  %.2e  %.2e  %.4e  %.2e  %.2e  %.2e  %.2e  %+.2e  %s\n', ...
                nit, evaluations.nfeval, iterate.feval,                     ...
                indicators.norm_glag, indicators.norm_z_s, vmax,            ...
                norm_d_x, norm_d_s, Delta_f, Delta_z, rho, it_type);
            if(~isempty(msg))
                disp(msg);
            end
        end
    else
        if (setting.show_errg)
            fprintf('%5d  %5d  %+.14e  %.2e  %.2e  %.4e %.2e                     %.2e  %+.2e\n',  ...
                nit, evaluations.nfeval, iterate.feval,                     ...
                indicators.norm_glag, indicators.norm_z_s, vmax, errg,      ...
                Delta_f, Delta_z);
             
        else
            fprintf('%5d  %5d  %+.14e  %.2e  %.2e  %.4e                      %.2e  %+.2e\n',      ...
                nit, evaluations.nfeval, iterate.feval,                     ...
                indicators.norm_glag, indicators.norm_z_s, vmax,            ...
                Delta_f, Delta_z);
        end
    end
end

end % end of deft_funnel_printout
