function deft_funnel_print_head_info( setting )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Prints the iteration info head.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
if ( setting.verbose >= 1 )
    
   fprintf( '\n  it    nfeval      fvalue            ' )

   if ( setting.show_errg )
      fprintf('optim   const_viol   vmax      errg    ||d_x||   ||d_s||')
   else
      fprintf('optim   const_viol   vmax      ||d_x||   ||d_s||')
   end
   
   fprintf('    Delta      rho      it_type\n')
   fprintf('\n')

end

end % end of deft_funnel_print_head_info