function [best_sol, best_feval, best_indicators, best_iterate, total_eval,  ...
    nb_local_searches, fL] = run_deft_funnel_multistart_single_bb_test_prob(nprob)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Runs DEFT-FUNNEL without multistart (global optimization) on a single 
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

[best_sol, best_feval, best_indicators, best_iterate, total_eval,           ...
    nb_local_searches, fL] = deft_funnel_multistart(f, c, h, dev_f, dev_h,  ...
    n, nb_cons_c, nb_cons_h, 'lsbounds', ls, 'usbounds', us, 'lxbounds',    ...
    lx, 'uxbounds', ux);

end
