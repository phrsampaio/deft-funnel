function c = problem_SR7_cons(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem Speed Reducer Design in "A Collection of Test Problems for 
% Constrained Global Optimization Algorithms", Springer-Verlag, Berlin,
% 1990, by C. A. Floudas and P. M. Pardalos.
%
% Desc: 
%     - Number of variables: 7
%     - Number of constraints (not bounds): 11 inequalities
%     - Objective function: non-linear
%     - Constraints: non-linear
%
% Lower and upper bounds for the constraint(s):
% lc = (-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf)
% uc = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
%
% Lower and upper bounds for the decision variables x:
% lx = (2.6, 0.7, 17, 7.3, 7.3, 2.9, 5)
% ux = (3.6, 0.8, 28, 8.3, 8.3, 3.9, 5.5)
%
% Global optimal sol (best known): 
% x* = (3.5, 0.7, 17, 7.3, 7.71, 3.35, 5.287)
% f(x*) = 2994.47
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A1 = ((745*x(4)/(x(2)*x(3)))^2 + (16.91*1000000))^(0.5);
B1 = 0.1*x(6)^3;
A2 = ((745*x(5)/(x(2)*x(3)))^2 + (157.5*1000000))^(0.5);
B2 = 0.1*x(7)^3;

c(1) = 27 - x(1)*(x(2)^2)*x(3);
c(2) = 397.5 - x(1)*(x(2)^2)*(x(3)^2);
c(3) = 1.93 - (x(2)*(x(6)^4)*x(3))/(x(4)^3);
c(4) = 1.93 - (x(2)*(x(7)^4)*x(3))/(x(5)^3);
c(5) = A1/B1 - 1100;
c(6) = A2/B2 - 850;
c(7) = x(2)*x(3) - 40;
c(8) = 5 - (x(1)/x(2));
c(9) = x(1)/x(2) - 12;
c(10) = 1.9 + 1.5*x(6) - x(4);
c(11) = 1.9 + 1.1*x(7) - x(5);

end