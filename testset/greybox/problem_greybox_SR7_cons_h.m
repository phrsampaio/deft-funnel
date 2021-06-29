function h = problem_greybox_SR7_cons_h(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem Speed Reducer Design in "A Collection of Test Problems for 
% Constrained Global Optimization Algorithms", Springer-Verlag, Berlin,
% 1990, by C. A. Floudas and P. M. Pardalos.
%
% Desc: 
%     - Number of variables: 7
%     - Number of constraints (not bounds): 11 inequalities:
%                                             - 9 black boxes (nonlinear)
%                                             - 2 white boxes (linear)
%     - Objective function: non-linear (white box)
%     - Constraints: 2 linear and 9 non-linear
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

h(1) = 1.9 + 1.5*x(6) - x(4);
h(2) = 1.9 + 1.1*x(7) - x(5);

end
