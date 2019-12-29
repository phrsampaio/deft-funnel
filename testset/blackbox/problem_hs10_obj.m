function f = problem_hs10_obj(x) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 10 from Hock and Schittkowski collection
%
% Desc: 
%     - Number of variables: 2
%     - Number of constraints (not bounds): 1 inequality
%     - Objective function: linear
%     - Constraints: non-linear
%
% Lower and upper bounds for the constraint(s):
% lc = 0
% uc = Inf
%
% Lower and upper bounds for the decision variables x:
% lx = (-Inf, -Inf)
% ux = (Inf, Inf)

% initial guess: x0 = (-10,10); f(x0) = -20 (not feasible)
% optimal sol:   x* = (0, 1);  f(x*) = -1
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = x(1)-x(2);

end