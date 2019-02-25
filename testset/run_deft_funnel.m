%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Runs DEFT-FUNNEL without multistart (local optimization) on a small 
% collection of test problems defined in:
% 1. 'deft_funnel_problem_init.m' (entry parameters),
% 2. 'deft_funnel_problem_obj.m' (objective function),
% 3. 'deft_funnel_problem_cons.m' (constraint function(s).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

for nprob = 1:10

   % Initialize test problem
   [ x0, nbcons, ls, us, lx, ux ] = deft_funnel_problem_init( nprob );

   disp( ' ' )
   disp( [' Running test problem #', int2str(nprob)] )
   
   [ x, fx, norm_c_s, mu, nfeval ] =                                     ...
       deft_funnel( @(x)deft_funnel_problem_obj( x, nprob ),             ...
       @(x)deft_funnel_problem_cons( x, nprob ), x0, nbcons,             ...
      'lsbounds', ls, 'usbounds', us, 'lxbounds', lx, 'uxbounds', ux )

end
