function f = problem_G6_obj(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem G6 in  "Evolutionary algorithms for constrained parameter 
% optimization problems.", Evolutionary Computation, 4(1):1–32, 1996,
% by Michalewicz, Z., Schoenauer, M.
%
% Desc: 
%     - Number of variables: 2
%     - Number of constraints (not bounds): 2
%     - Objective function: non-linear
%     - Constraints: non-linear
%
% Lower and upper bounds for the constraint(s):
% lc = (0, 0)
% uc = (Inf, Inf)
%
% Lower and upper bounds for the decision variables x:
% lx = (13, 0)
% ux = (100, 100)
%
% Initial guess:      x0 = (20.1, 5.84) (not feasible)
% Global optimal sol: x* = (14.095, 0.84296); f(x*) = -6961.81381
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = (x(1)-10)^3 + (x(2)-20)^3;

end