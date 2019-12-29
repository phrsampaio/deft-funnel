function output = problem_greybox_hesse_dev_h(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: "A Heuristic Search Procedure for Estimating a Global Solution of 
% Nonconvex Programming Problems", Operations Research, 21(6):1267--1280, 1973,
% by Rick Hesse.
%
% Category:
%     - Objective function: concave quadratic (white box)
%     - Constraints: 3 linear (black boxes) and 2 quadratic (white boxes)
%     - Type: 3 separable problems combined
%     - 18 local minima and 1 global minimum
%
% Initial guess: not given
% Global optimal sol: x* = (5,1,5,0,5,10); f(x*) = -310
% lb = (0, 0, 1, 0, 1, 0)
% ub = (Inf, Inf, 5, 6, 5, 10)
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Derivatives of the white-box constraints
Jh = [0 0 2*(x(3)-3) 1 0 0; 0 0 0 0 2*(x(5)-3) 1];
Hh1 = [zeros(1,6); zeros(1,6); 0 0 2 0 0 0; zeros(1,6); zeros(1,6); zeros(1,6)];
Hh2 = [zeros(1,6); zeros(1,6); zeros(1,6); zeros(1,6); 0 0 0 0 2 0; zeros(1,6)];
Hhlist = [Hh1 Hh2];
output = {Jh, Hhlist};

end
