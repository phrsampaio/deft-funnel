function [x, fx, mu, indicators, evaluations, iterate, exit_algo] =         ...
    run_deft_funnel_single_bb_test_prob(nprob)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Runs DEFT-FUNNEL without multistart (local optimization) on a single 
% test problem defined in:
% 1. 'deft_funnel_problem_init.m' (entry parameters),
% 2. 'deft_funnel_problem_obj.m' (objective function),
% 3. 'deft_funnel_problem_cons_c.m' (constraint function(s)).
% The test problem to be solved is defined by the input 'nprob'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x0, n, nb_cons_c, nb_cons_h, ls, us, lx, ux] = deft_funnel_problem_init(nprob);

f = @(x)deft_funnel_problem_obj(x, nprob);
c = @(x)deft_funnel_problem_cons_c(x, nprob);
h = [];
dev_f = [];
dev_h = [];

[x, fx, mu, indicators, evaluations, iterate, exit_algo] =                  ...
    deft_funnel(f, c, h, dev_f, dev_h, x0, nb_cons_c, nb_cons_h,            ...
    'lsbounds', ls, 'usbounds', us, 'lxbounds', lx, 'uxbounds', ux)

end
