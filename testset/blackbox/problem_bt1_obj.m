function f = problem_bt1_obj(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: problem 1 from P.T. Boggs and J.W. Tolle,
% "A strategy for global convergence in a sequential 
% quadratic programming algorithm", SINUM 26(3), pp. 600-623, 1989.
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
% Initial guess: x0 = (0.08,0.06)
% Optimal sol:   x* = (1,0)
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = 10*x(1)^2+10*x(2)^2-x(1)-10;

end
