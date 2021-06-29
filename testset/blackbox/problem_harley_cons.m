function c = problem_harley_cons(x)

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

% Translation from the source
% x -> x(1)
% y -> x(2)
% A -> x(3)
% B -> x(4)
% Cx -> x(5)
% Cy -> x(6)
% Px -> x(7)
% Py -> x(8)
% p -> x(9)

c(1) = x(7) + x(8) - x(3) - x(4);
c(2) = x(1) - x(7) - x(5);
c(3) = x(2) - x(8) - x(6);
c(4) = x(9)*x(7) + 2*x(5) - 2.5*x(1);
c(5) = x(9)*x(8) + 2*x(6) - 1.5*x(2);
c(6) = x(9)*x(7) + x(9)*x(8) - 3*x(3) - x(4);

end
