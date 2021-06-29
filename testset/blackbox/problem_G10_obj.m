function f = problem_G10_obj(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem G10 in "Evolutionary algorithms for constrained parameter 
% optimization problems.", Evolutionary Computation, 4(1):1–32, 1996,
% by Michalewicz, Z., Schoenauer, M.
% Also the problem 106 (heat exchanger design) from Hock and Schittkowski 
% collection.
%
% Desc: 
%     - Number of variables: 8
%     - Number of constraints (not bounds): 6 inequalities
%     - Objective function: linear
%     - Constraints: 3 linear and 3 non-linear
%     - All constraints are active at the optimum.
%
% Lower and upper bounds for the constraint(s):
% lc = (0, 0, 0, 0, 0, 0)
% uc = (Inf, Inf, Inf, Inf, Inf, Inf)
%
% Lower and upper bounds for the decision variables x:
% lx = (100, 1000, 1000, 10, 10, 10, 10, 10)
% ux = (10000, 10000, 10000, 1000, 1000, 1000, 1000, 1000)
%
% Initial guess: 
% x0 = (5000, 5000, 5000, 200, 350, 150, 225, 425)
% f(x0) = 15000 (not feasible)
%
% Global optimal sol: 
% x* = (579.3167, 1359.943, 5110.071, 182.0174, 295.5985, 217.9799, 286.4162, 395.5979)
% f(x*) = 7049.330923.
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = x(1) + x(2) + x(3);

end