function [ x, fx, mu, indicators, evaluations ] = run_deft_funnel_single_prob( nprob )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Runs DEFT-FUNNEL without multistart (local optimization) on a single 
% test problem defined in:
% 1. 'deft_funnel_problem_init.m' (entry parameters),
% 2. 'deft_funnel_problem_obj.m' (objective function),
% 3. 'deft_funnel_problem_cons.m' (constraint function(s).
% The test problem to be solved is defined by the input 'nprob'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ x0, nbcons, ls, us, lx, ux ] = deft_funnel_problem_init( nprob );

[x, fx, mu, indicators, evaluations ] =                                     ...
   deft_funnel( @(x)deft_funnel_problem_obj( x, nprob ),                    ...
   @(x)deft_funnel_problem_cons( x, nprob ), x0, nbcons,                    ...
  'lsbounds', ls, 'usbounds', us, 'lxbounds', lx, 'uxbounds', ux );

end