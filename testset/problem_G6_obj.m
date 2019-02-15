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
% Initial guess:      x0 = (20.1, 5.84)
% Global optimal sol: x* = (14.095, 0.84296); f(x*) = -6961.81381
% lb = (13, 0)
% ub = (100, 100)
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = (x(1)-10)^3 + (x(2)-20)^3;

end