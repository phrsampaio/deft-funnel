function c = problem_handbook_quadcons_pb3_cons( x )

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

c(1) = (x(3)-3)^2 + x(4);
c(2) = (x(5)-3)^2 + x(6);
c(3) = x(1) - 3*x(2);
c(4) = -x(1) + x(2);
c(5) = x(1) + x(2);

end