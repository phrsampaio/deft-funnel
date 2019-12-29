function c = problem_G9_cons(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem G9 in  "Evolutionary algorithms for constrained parameter 
% optimization problems.", Evolutionary Computation, 4(1):1–32, 1996,
% by Michalewicz, Z., Schoenauer, M.
%
% Desc: 
%     - Number of variables: 7
%     - Number of constraints (not bounds): 4 inequalities
%     - Objective function: non-linear
%     - Constraints: non-linear
%     - Type: First and last constraints active in the global optimum
%
% Lower and upper bounds for the constraint(s):
% lc = (0, 0, 0, 0)
% uc = (Inf, Inf, Inf, Inf)
%
% Lower and upper bounds for the decision variables x:
% lx = -10, i=1..7
% ux = 10, i=1..7
%
% Global optimal sol:   
% x*    = (2.330499, 1.951372, -0.4775414, 4.365726, -0.6244870, 1.038131, 1.594227)
% f(x*) = 680.6300573
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c(1) = 127 - 2*x(1)^2 - 3*x(2)^4 - x(3) - 4*x(4)^2 - 5*x(5);
c(1) = c(1)/127;
c(2) = 282 - 7*x(1) - 3*x(2) - 10*x(3)^2 - x(4) + x(5);
c(2) = c(2)/282;
c(3) = 196 - 23*x(1) - x(2)^2 - 6*x(6)^2 + 8*x(7);
c(3) = c(3)/196;
c(4) = -4*x(1)^2 - x(2)^2 + 3*x(1)*x(2) - 2*x(3)^2 - 5*x(6) + 11*x(7);

end