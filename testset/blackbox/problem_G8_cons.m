function c = problem_G8_cons(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem G8 in  "Evolutionary algorithms for constrained parameter 
% optimization problems.", Evolutionary Computation, 4(1):1–32, 1996,
% by Michalewicz, Z., Schoenauer, M.
%
% Desc: 
%     - Number of variables: 2
%     - Number of constraints (not bounds): 2 inequalities
%     - Objective function: non-linear
%     - Constraints: non-linear
%     - Type: many local minima with the lowest peaks located along the x axis
%     - Maximization problem turned into mim. problem
%
% Lower and upper bounds for the constraint(s):
% lc = (-Inf, -Inf)
% uc = (0, 0)
%
% Lower and upper bounds for the decision variables x:
% lx = (0, 0)
% ux = (10, 10)
%
% Global optimal sol: x* = (1.2279713, 4.2453733); f(x*) = -0.095825
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c(1) = x(1)^2 - x(2) + 1;
c(2) = 1 - x(1) + (x(2) - 4)^2;

end