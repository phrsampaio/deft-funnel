function f = problem_hs13_obj( x ) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem 13 from Hock and Schittkowski collection
% Initial guess: x0 = (-2,-2); f(x0) = 20 (not feasible)
% Optimal sol:   x* = (1, 0);  f(x*) = 1
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = (x(1)-2)^2 + x(2)^2;

end