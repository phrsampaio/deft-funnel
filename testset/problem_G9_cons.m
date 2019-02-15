function c = problem_G9_cons(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem G9 in  "Evolutionary algorithms for constrained parameter 
% optimization problems.", Evolutionary Computation, 4(1):1–32, 1996,
% by Michalewicz, Z., Schoenauer, M.
%
% Desc: 
%     - Number of variables: 7
%     - Number of constraints (not bounds): 4
%     - Objective function: non-linear
%     - Constraints: non-linear
%     - Type: First and last constraints active in the global optimum
%
% Global optimal sol:   
% x*    = (2.330499, 1.951372, -0.4775414, 4.365726, -0.6244870, 1.038131, 1.594227)
% f(x*) = 680.6300573
% lb    = -10, i=1..7
% ub    = 10, i=1..7
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v1 = 2*x(1)^2;
v2 = x(2)^2;
c(1) = v1 + 3*v2^2 + x(3) + 4*x(4)^2 + 5*x(5) - 127;
c(2) = 7*x(1) + 3*x(2) + 10*x(3)^2 + x(4) - x(5) - 282;
c(3) = 23*x(1) + v2 + 6*x(6)^2 - 8*x(7) - 196;
c(4) = 2*v1 + v2 - 3*x(1)*x(2) + 2*x(3)^2 + 5*x(6) - 11*x(7);

end