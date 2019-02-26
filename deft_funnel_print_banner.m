function deft_funnel_print_banner( setting )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Prints the algorithm's banner.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( setting.verbose >= 1 )

   fprintf( '\n' )
   fprintf( '  ********************************************************' )
   fprintf( '***************************************\n' )
   fprintf( '  *                                                       ' )
   fprintf( '                                      *\n' )
   fprintf( '  *                        DEFT-FUNNEL: Derivative-Free ' )
   fprintf( 'Trust FUNNEL                            *\n' )
   fprintf( '  *                                                       ' )
   fprintf( '                                      *\n' )
   fprintf( '  *                          Copyright (c) 2019 ' )
   fprintf( 'Phillipe R. Sampaio                             *\n')
   fprintf( '  *                                                       ' )
   fprintf( '                                      *\n')
   fprintf( '  ********************************************************' )
   fprintf( '***************************************\n' )

end

end % end of deft_funnel_print_banner