function f = problem_handbook_quadcons_pb3_obj( x )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Test problem 3 in "Quadratically Constrained Problems" section of
% "Handbook of test problems in local and Global Optimization", pp.
% 24-25, by Carl A. Schweiger, Christodoulos A. Floudas, Claire Adjiman, 
% Clifford A. Meyer, John L. Klepeis, Panos M. Pardalos, Stephen T. Harding, 
% William R. Esposito, and Zeynep H. GümüsSINUM, 1999.
%
% Category:
%     - Objective function: concave quadratic
%     - Constraints: linear and quadratic
%     - Type: 3 separable problems combined
%     - 18 local minima and 1 global minimum
%
% Initial guess: not given
% Global optimal sol: x* = (5,1,5,0,5,10); f(x*) = -310
% lb = (0, 0, 1, 0, 1, 0)
% ub = (Inf, Inf, 5, 6, 5, 10)
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = -25*(x(1)-2)^2 - (x(2)-2)^2 - (x(3)-1)^2 - (x(4)-4)^2 - ...
     (x(5)-1)^2 - (x(6)-4)^2;

end
