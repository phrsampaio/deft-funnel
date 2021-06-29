function c = problem_greybox_WB4_cons_h(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Welded Beam Design Problem in "An efficient constraint handling 
% method for genetic algorithms", Computer Methods in Applied Mechanics and 
% Engineering 186, 311–338, 2002, by Deb, K.
%
% Desc: 
%     - Number of variables: 4
%     - Number of constraints (not bounds): 6 inequalities
%     - Objective function: non-linear
%     - Constraints: non-linear
%     - Type: engineering design problem
%
% Lower and upper bounds for the constraint(s):
% lc = (-Inf, -Inf, -Inf, -Inf, -Inf, -Inf)
% uc = (0, 0, 0, 0, 0, 0)
%
% Lower and upper bounds for the decision variables x:
% lx = (0.125, 0.1, 0.1, 0.1)
% ux = (10, 10, 10, 10)
%
% Global optimal sol: 
% x*    = (0.20564426101885, 3.47257874213172, 9.03662391018928, 0.20572963979791)
% f(x*) = 1.7250
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set constants
xmax = 10;

% Define the constraints
c(1) = (x(1) - x(4))/xmax;
c(2) = (0.10471*x(1)^2 + 0.04811*x(3)*x(4)*(14.0 + x(2)) - 5.0)/5.0;

end