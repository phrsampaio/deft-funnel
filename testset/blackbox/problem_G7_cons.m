function c = problem_G7_cons(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem G7 in  "Evolutionary algorithms for constrained parameter 
% optimization problems.", Evolutionary Computation, 4(1):1–32, 1996,
% by Michalewicz, Z., Schoenauer, M.
%
% Desc: 
%     - Number of variables: 10
%     - Number of constraints (not bounds): 8 inequalties
%     - Objective function: quadratic
%     - Constraints: 3 linear and 5 non-linear
%     - All constraints are active at the optimum except the last two
%
% Lower and upper bounds for the constraint(s):
% lc = (0, 0, 0, 0, 0, 0, 0, 0)
% uc = (Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf)
%
% Lower and upper bounds for the decision variables x:
% lx = (-10, -10, -10, -10, -10, -10, -10, -10, -10, -10)
% ux = (10, 10, 10, 10, 10, 10, 10, 10, 10, 10)
%
% Global optimal sol: 
% x* = (2.171996, 2.363683, 8.773926, 5.095984, 0.9906548,
%       1.430574, 1.321644, 9.828726, 8.280092, 8.375927)
% f(x*) = 24.3062091
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c(1) = 105 - 4*x(1) - 5*x(2) + 3*x(7) - 9*x(8);
c(1) = c(1)/105;
c(2) = -3*(x(1)-2)^2 - 4*(x(2)-3)^2 - 2*x(3)^2 + 7*x(4) + 120;
c(2) = c(2)/1258;
c(3) = -10*x(1) + 8*x(2) + 17*x(7) - 2*x(8);
c(3) = c(3)/370;
c(4) = -x(1)^2 - 2*(x(2)-2)^2 + 2*x(1)*x(2) - 14*x(5) + 6*x(6);
c(4) = c(4)/788;
c(5) = 8*x(1) - 2*x(2) - 5*x(9) + 2*x(10) + 12;
c(5) = c(5)/158;
c(6) = -5*x(1)^2 - 8*x(2) - (x(3)-6)^2 + 2*x(4) + 40;
c(6) = c(6)/816;
c(7) = 3*x(1) - 6*x(2) - 12*(x(9)-8)^2 + 7*x(10);
c(7) = c(7)/4048;
c(8) = -0.5*(x(1)-8)^2 - 2*(x(2)-4)^2 - 3*x(5)^2 + x(6) + 30;
c(8) = c(8)/834;

end