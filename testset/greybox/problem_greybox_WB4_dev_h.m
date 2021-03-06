function output = problem_greybox_WB4_dev_h(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Welded Beam Design Problem in "An efficient constraint handling 
% method for genetic algorithms", Computer Methods in Applied Mechanics and 
% Engineering 186, 311�338, 2002, by Deb, K.
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

% Derivatives of the white-box constraints
Jh = [1/xmax 0 0 -1/xmax; ...
      0.10471*2*x(1)/5.0 0.04811*x(3)*x(4)/5.0 ...
      0.04811*x(4)*(14.0 + x(2))/5.0 0.04811*x(3)*(14.0 + x(2))/5.0];
Hh1 = [zeros(1,4); zeros(1,4); zeros(1,4); zeros(1,4)];

dx1dx = [0.10471*2/5.0 0 0 0];
dx2dx = [0 0 0 0];
dx3dx = [0 0.04811*x(4)/5.0 0 0.04811*(14.0 + x(2))/5.0];
dx4dx = [0 0.04811*x(2)/5.0 0.04811*x(3)*(14.0 + x(2))/5.0 0];
    
Hh2 = [dx1dx; dx2dx; dx3dx; dx4dx];
Hhlist = [Hh1 Hh2];
output = {Jh, Hhlist};

end