function deft_funnel_printout( nit, evaluations, iterate, setting, ...
    indicators, vmax, norm_d_x, norm_d_s, Delta, rho, it_type, msg, errg )
      
      if ( setting.verbose >= 1 )
         if ( nit > 0 )
             if ( setting.show_errg )
                fprintf( '%5d  %5d  %+.14e  %.2e  %.2e  %.4e %.2e %.2e  %.2e  %.2e  %+.2e  %s\n', ...
                    nit, evaluations.nfeval, iterate.feval, indicators.norm_glag, indicators.norm_c_s, vmax,          ...
                    errg, norm_d_x, norm_d_s, Delta, rho, it_type );
                if( ~isempty( msg ) )
                    disp( msg );
                end
             else
                fprintf( '%5d  %5d  %+.14e  %.2e  %.2e  %.4e  %.2e  %.2e  %.2e  %+.2e  %s\n', ...
                    nit, evaluations.nfeval, iterate.feval, indicators.norm_glag, indicators.norm_c_s, vmax,          ...
                    norm_d_x, norm_d_s, Delta, rho, it_type );
                if( ~isempty( msg ) )
                    disp( msg );
                end
             end
         else
             if ( setting.show_errg )
                fprintf( '%5d  %5d  %+.14e  %.2e  %.2e  %.4e %.2e                     %.2e\n',  ...
                    nit, evaluations.nfeval, iterate.feval, indicators.norm_glag, indicators.norm_c_s, ...
                    vmax, errg, Delta );
                 
             else
                fprintf( '%5d  %5d  %+.14e  %.2e  %.2e  %.4e                      %.2e\n',      ...
                    nit, evaluations.nfeval, iterate.feval, indicators.norm_glag, indicators.norm_c_s, vmax, Delta );
             end
         end
      end
end
