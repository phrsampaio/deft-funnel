function f = problem_hs12_obj( x ) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem 12 from Hock and Schittkowski collection
% Initial guess: x0 = (0,0); f(x0) = 0 (feasible)
% Optimal sol:   x* = (2, 3);  f(x*) = -30
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = 0.5*x(1)^2 + x(2)^2 - x(1)*x(2) - 7*x(1) - 7*x(2);

end