function f = problem_gomez_pb3_obj(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Test problem 3 in "The tunnelling method for solving the 
% constrained global optimization problem with several non-connected 
% feasible regions",
% by Susana Gómez and A. V. Levy, In book: Numerical Analysis, pp.34-47, 1982.
%
% Desc: 
%     - Number of variables: 2
%     - Number of constraints (not bounds): 1 inequality
%     - Objective function: non-convex 
%     - Constraint: non-convex 
%     - Type: disconnected feasible region
%     - 1 global minimum
%
% Lower and upper bounds for the constraint(s):
% lc = -Inf
% uc = 0
%
% Lower and upper bounds for the decision variables x:
% lx = (-1, -1)
% ux = (1, 1)
%  
% Initial guess: not given
% Global optimal sol: x* = (0.109, -0.623); f(x*) = -0.9711
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = (4 -2.1*x(1)^2 + (x(1)^4)/3)*(x(1)^2) + x(1)*x(2) + (-4 + 4*x(2)^2)*(x(2)^2);

end
