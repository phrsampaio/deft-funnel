function c = problem_PrP_cons(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Pressure Vessel Design problem in "Constraint-handling in genetic 
% algorithms through the use of dominance-based tournament selection", 
% Advanced Engineering Informatics, 16, 193–203, 2002, by 
% Coello Coello, C. A. and Mezura-Montes, E..
%
% Desc: 
%     - Number of variables: 4
%     - Number of constraints (not bounds): 3
%     - Objective function: nonlinear
%     - Constraint: 1 linear and 2 nonlinear
%     - Type: engineering design problem
%  
% Best known value: f(x*) = 5804.45
% lb = (0, 0, 0, 0)
% ub = (1, 1, 50, 240)
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c(1) = -x(1) + 0.0193*x(3);
c(2) = -x(2) + 0.00954*x(3);
c(3) = -pi*x(3)^2*x(4) - (4/3)*pi*x(3)^3 + 1296000;

end