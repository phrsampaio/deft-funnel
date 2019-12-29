function f = problem_hs7_obj(x) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem 7 from Hock and Schittkowski collection
%
% Desc: 
%     - Number of variables: 2
%     - Number of constraints (not bounds): 1 equality
%     - Objective function: non-linear
%     - Constraints: non-linear
%
% Lower and upper bounds for the constraint(s):
% lc = 0
% uc = 0
%
% Lower and upper bounds for the decision variables x:
% lx = (-Inf, -Inf)
% ux = (Inf, Inf)
%
% Initial guess: x0 = (2,2); f(x0) = log(5)-2 (not feasible)
% Optimal sol:   x* = (0,sqrt(3));  f(x*) = -sqrt(3)
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = log(1+x(1)^2)-x(2);

end