function c = problem_gomez_pb3_cons( x )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Test problem 3 in "The tunnelling method for solving the 
% constrained global optimization problem with several non-connected 
% feasible regions",
% by Susana G�mez and A. V. Levy, In book: Numerical Analysis, pp.34-47, 1982.
%
% Desc: 
%     - Number of variables: 2
%     - Number of constraints (not bounds): 1
%     - Objective function: non-convex 
%     - Constraint: non-convex 
%     - Type: disconnected feasible region
%     - 1 global minimum
%  
% Initial guess: not given
% Global optimal sol: x* = (0.109, -0.623); f(x*) = -0.9711
% lb = (-1, -1)
% ub = (1, 1)
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = -sin(4*pi*x(1)) + 2*sin(2*pi*x(2))^2;

end