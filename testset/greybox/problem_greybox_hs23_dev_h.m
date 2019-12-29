function output = problem_greybox_hs23_dev_h(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem 23 from Hock and Schittkowski collection
%
% Desc: 
%     - Number of variables: 2
%     - Number of constraints (not bounds): 5 inequalities (2 black boxes and
%                                           3 white boxes)
%     - Objective function: non-linear (black-box)
%     - Constraints: 1 linear and 4 non-linear
%
% Lower and upper bounds for the constraint(s):
% lc = (0, 0, 0, 0, 0)
% uc = (Inf, Inf, Inf, Inf, Inf)
%
% Lower and upper bounds for the decision variables x:
% lx = (-50, -50)
% ux = (50, 50)
%
% Initial guess: x0 = (3,1); f(x0) = 10 (not feasible)
% Optimal sol:   x* = (1,1); f(x*) = 2
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Derivatives of the white-box constraints
Jh = [18*x(1) 2*x(2); 2*x(1) -1; -1 2*x(2)];
Hh1 = [18 0; 0 2];
Hh2 = [2 0; 0 0];
Hh3 = [0 0; 0 2];
Hhlist = [Hh1 Hh2 Hh3];
output = {Jh, Hhlist};

end
