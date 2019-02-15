clear all

for nprob = 1:10

   % Initialize test problem
   [ x0, nbcons, ls, us, lx, ux ] = deft_funnel_init( nprob );

   disp( ' ' )
   disp( [' Running test problem #', int2str(nprob)] )
   
   [ x, fx, norm_c_s, mu, nfeval ] =                                     ...
       deft_funnel( @(x)deft_funnel_problem_obj( x, nprob ),             ...
       @(x)deft_funnel_problem_cons( x, nprob ), x0, nbcons,             ...
      'lsbounds', ls, 'usbounds', us, 'lxbounds', lx, 'uxbounds', ux )

end
