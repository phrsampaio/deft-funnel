function deft_funnel_print_head_info( setting )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Prints the iteration info head.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
fprintf( '\n')
fprintf( '  it    nfeval      fvalue            ' )
   
if ( setting.verbose >= 1 )

   if ( setting.show_errg )
      fprintf('optim   const_viol   vmax      errg    ||d_x||   ||d_s||')
   else
      fprintf('optim   const_viol   vmax      ||d_x||   ||d_s||')
   end
   
   fprintf('    Delta      rho      it_type\n')
   fprintf('\n')

end

end % end of deft_funnel_print_head_info