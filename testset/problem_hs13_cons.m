function c = problem_hs13_cons( x ) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem 13 from Hock and Schittkowski collection
% Initial guess: x0 = (-2,-2); f(x0) = 20 (not feasible)
% Optimal sol:   x* = (1, 0);  f(x*) = 1
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = (1-x(1))^3 - x(2);

end
