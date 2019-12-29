function c = problem_hs12_cons(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem 12 from Hock and Schittkowski collection
%
% Desc: 
%     - Number of variables: 2
%     - Number of constraints (not bounds): 1 inequality
%     - Objective function: non-linear
%     - Constraints: non-linear
%
% Lower and upper bounds for the constraint(s):
% lc = 0
% uc = Inf
%
% Lower and upper bounds for the decision variables x:
% lx = (-Inf, -Inf)
% ux = (Inf, Inf)
%
% Initial guess: x0 = (0,0); f(x0) = 0 (feasible)
% Optimal sol:   x* = (2, 3);  f(x*) = -30
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 25 - 4*x(1)^2 - x(2)^2;

end
