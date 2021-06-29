function c = problem_G11_cons(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem G11 in  "Evolutionary algorithms for constrained parameter 
% optimization problems.", Evolutionary Computation, 4(1):1–32, 1996,
% by Michalewicz, Z., Schoenauer, M.
%
% Desc: 
%     - Number of variables: 2
%     - Number of constraints (not bounds): 1 equality
%     - Objective function: non-linear
%     - Constraint: non-linear
%
% Lower and upper bounds for the constraint(s):
% lc = 0
% uc = 0
%
% Lower and upper bounds for the decision variables x:
% lx = -1, i=1..2
% ux = 1, i=1..2
%
% Global optimal sol: x* = (+- 0.70711, 0.5); f(x*) = 0.75000455
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = x(2)-x(1)^2;

end