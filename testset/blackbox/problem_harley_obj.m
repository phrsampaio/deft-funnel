function f = problem_harley_obj(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Harley Pooling Problem (5.2.2) in "Handbook of Test Problems in 
% Local and Global, Nonconvex optimization and its applications", 
% Springer-Verlag US, 1999, by C. A. Floudas et al.
%
% Desc: 
%     - Number of variables: 9
%     - Number of constraints (not bounds): 3 linear equalities
%                                           1 nonlinear equality
%                                           2 nonconvex inequalties
%     - Objective function: non-linear
%     - Type: pooling problem
%
% Lower and upper bounds for the constraint(s):
% lc = (0, 0, 0, -Inf, -Inf, 0)
% uc = (0, 0, 0, 0, 0, 0)
%
% Lower and upper bounds for the decision variables x:
% lx = (0, 0, 0, 0, 0, 0, 0, 0, 0)
% ux = (600, 200, 500, 500, 500, 500, 500, 500, 500)
%
% Global optimal sol: 
% x*    = (600, 0, 300, 0, 300, 0, 300, 0, 3)
% f(x*) = -600
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = -9*x(1) - 15*x(2) + 6*x(3) + 16*x(4) + 10*(x(5) + x(6));

end
