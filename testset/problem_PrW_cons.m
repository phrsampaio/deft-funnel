function c = problem_PrW_cons(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Welded Beam Design Problem in "An efficient constraint handling 
% method for genetic algorithms", Computer Methods in Applied Mechanics and 
% Engineering 186, 311–338, 2002, by Deb, K.
%
% Desc: 
%     - Number of variables: 4
%     - Number of constraints (not bounds): 4
%     - Objective function: non-linear
%     - Constraints: non-linear
%     - Type: engineering design problem
%
% Global optimal sol obtained by the FSA algorithm: 
% x*    = (0.20564426101885, 3.47257874213172, 9.03662391018928, 0.20572963979791)
% f(x*) = 1.728226 (new:  1.7250022 by FSA)
% lb    = (0.125, 0.1, 0.1, 0.1)
% ub    = (10, 10, 10, 10)
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set constants
P = 6000; 
L = 14; 
E = 30e+6; 
G = 12e+6;

tmax = 13600; 
smax = 30000; 
dmax = 0.25;

M  = P * (L + x(2)/2); 
R  = sqrt(0.25 * (x(2)^2 + (x(1) + x(3))^2));
J  = sqrt(2)* x(1) * x(2) * (x(2)^2/12 + 0.25 * (x(1) + x(3))^2);
Pc = (4.013*E / (6*L^2)) * x(3) *x(4)^3 * (1 - 0.25 * x(3) * sqrt(E/G)/L);
t1 = P/(sqrt(2)*x(1)*x(2)); 
t2 = M*R/J;
t  = sqrt(t1^2 + t1*t2*x(2)/R + t2^2);
s  = 6*P*L / (x(4) * x(3)^2);
d  = 4*P*L^3 / (E * x(4) * x(3)^3);

% Define the constraints
c(1) = t - tmax;
c(2) = s - smax;
c(3) = x(1) - x(4);
c(4) = 0.10471*x(1)^2 + 0.04811*x(3)*x(4)*(14.0 + x(2)) - 5.0;
c(5) = d - dmax;
c(6) = P - Pc;

end