function f = problem_G8_obj(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem G8 in  "Evolutionary algorithms for constrained parameter 
% optimization problems.", Evolutionary Computation, 4(1):1–32, 1996,
% by Michalewicz, Z., Schoenauer, M.
%
% Desc: 
%     - Number of variables: 2
%     - Number of constraints (not bounds): 2
%     - Objective function: non-linear
%     - Constraints: non-linear
%     - Type: many local minima with the lowest peaks located along the x axis
%     - Maximization problem turned into mim. problem
%
% Global optimal sol: x* = (1.2279713, 4.2453733); f(x*) = -0.095825
% lb = (0, 0)
% ub = (10, 10)
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = -( sin(2*pi*x(1))^3 * sin(2*pi*x(2)) ) / ( x(1)^3 * (x(1)+x(2)) );

end