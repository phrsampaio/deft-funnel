%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Runs DEFT-FUNNEL without multistart (local optimization) on a small 
% collection of black-box test problems defined in:
% 1. 'deft_funnel_problem_init.m' (entry parameters),
% 2. 'deft_funnel_problem_obj.m' (objective function),
% 3. 'deft_funnel_problem_cons_c.m' (constraint function(s)).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

for nprob = 1:9

    % Initialize test problem
    [x0, n, nb_cons_c, nb_cons_h, ls, us, lx, ux] = deft_funnel_problem_init(nprob);
    f = @(x)deft_funnel_problem_obj(x, nprob);
    c = @(x)deft_funnel_problem_cons_c(x, nprob);
    h = [];
    dev_f = [];
    dev_h = [];
   
    disp(' ')
    disp([' Running test problem #', int2str(nprob)])
   
    [x, fx, mu, indicators, evaluations, iterate, exit_algo] =              ...
        deft_funnel(f, c, h, dev_f, dev_h, x0, nb_cons_c, nb_cons_h,        ...
        'lsbounds', ls, 'usbounds', us, 'lxbounds', lx, 'uxbounds', ux)

end
