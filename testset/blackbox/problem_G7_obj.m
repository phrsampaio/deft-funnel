function f = problem_G7_obj(x)

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

f = x(1)^2 + x(2)^2 + x(1)*x(2) - 14*x(1) - 16*x(2) + (x(3)-10)^2 +         ...
    4*(x(4)-5)^2 + (x(5)-3)^2 + 2*(x(6)-1)^2 + 5*x(7)^2 + 7*(x(8)-11)^2 +   ...
    2*(x(9)-10)^2 + (x(10)-7)^2 + 45;

end