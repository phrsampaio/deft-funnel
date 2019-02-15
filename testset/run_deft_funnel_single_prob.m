function [ x, fx, mu, indicators, evaluations ] = run_deft_funnel_single_prob( nprob )

[ x0, nbcons, ls, us, lx, ux ] = deft_funnel_init( nprob );

[x, fx, mu, indicators, evaluations ] =                                     ...
   deft_funnel( @(x)deft_funnel_problem_obj( x, nprob ),                    ...
   @(x)deft_funnel_problem_cons( x, nprob ), x0, nbcons,                    ...
  'lsbounds', ls, 'usbounds', us, 'lxbounds', lx, 'uxbounds', ux );

end