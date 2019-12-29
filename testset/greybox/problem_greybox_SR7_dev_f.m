function output = problem_greybox_SR7_dev_f(x)

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

A = 3.3333*x(3)^2 + 14.9334*x(3) - 43.0934;
B = x(6)^2 + x(7)^2;
C = x(6)^3 + x(7)^3;
D = x(4)*x(6)^2 + x(5)*x(7)^2;

% Derivative of the white-box objective function
gf = [(0.7854*A*x(2)^2 - 1.508*B)                                           ...
    0.7854*A*x(1)*2*x(2)                                                    ...
    (0.7854*x(1)*x(2)^2*(3.3333*2*x(3) + 14.9334))                          ...
    0.7854*x(6)^2                                                           ...
    0.7854*x(7)^2                                                           ...
    (-1.508*x(1)*2*x(6) + 7.477*3*x(6)^2 + 0.7854*x(4)*2*x(6))              ...
    (-1.508*x(1)*2*x(7) + 7.477*3*x(7)^2 + 0.7854*x(5)*2*x(7))]';
    
dx1dx = [0 0.7854*A*2*x(2) 0.7854*x(2)^2*(3.3333*2*x(3) + 14.9334)          ...
    0 0 -1.508*2*x(6) -1.508*2*x(7)];
    
dx2dx = [0.7854*A*2*x(2) 0.7854*A*x(1)*2 0 0 0 0 0];

dx3dx = [(0.7854*x(2)^2*(3.3333*2*x(3) + 14.9334))                          ...
    (0.7854*x(1)*2*x(2)*(3.3333*2*x(3) + 14.9334))                          ...
    (0.7854*x(1)*x(2)^2*3.3333*2)                                           ...
    0 0 0 0];
    
dx4dx = [0 0 0 0 0 0.7854*2*x(6) 0];

dx5dx = [0 0 0 0 0 0 0.7854*2*x(7)];

dx6dx = [-1.508*2*x(6) 0 0 0.7854*2*x(6) 0 (-1.508*x(1)*2 + 0.7854*x(4)*2) 0];

dx7dx = [-1.508*2*x(7) 0 0 0 0.7854*2*x(7) 0 (-1.508*x(1)*2 + 0.7854*x(5)*2)];
    
Hf = [dx1dx; dx2dx; dx3dx; dx4dx; dx5dx; dx6dx; dx7dx];
output = {gf, Hf};

end
