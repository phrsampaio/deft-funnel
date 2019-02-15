function c = problem_G3_cons(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem G3 in  "Evolutionary algorithms for constrained parameter 
% optimization problems.", Evolutionary Computation, 4(1):1–32, 1996,
% by Michalewicz, Z., Schoenauer, M.
%
% Desc: 
%     - Number of variables: n
%     - Number of constraints (not bounds): 1
%     - Objective function: non-linear
%     - Constraints: non-linear
%     - Global minimum in the boundary
%
% Global optimal sol: x* = (1/sqrt(n), ..., 1/sqrt(n)); f(x*) = -1
% lb = (0, 0)
% ub = (1, 1)
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = sum(x.^2)-1;

end