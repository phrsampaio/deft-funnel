function f = problem_G11_obj(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem G11 in  "Evolutionary algorithms for constrained parameter 
% optimization problems.", Evolutionary Computation, 4(1):1–32, 1996,
% by Michalewicz, Z., Schoenauer, M.
%
% Desc: 
%     - Number of variables: 2
%     - Number of inequalities (not bounds): 0
%     - Number of equalities (not bounds): 1
%     - Objective function: non-linear
%     - Constraint: non-linear
%
% Global optimal sol: x* = (+- 0.70711, 0.5); f(x*) = 0.75000455
% lb = -1, i=1..2
% ub = 1, i=1..2
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = x(1)^2 + (x(2)-1)^2;

end