function c = problem_hs7_cons( x ) % Problem 7 HS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem 7 from Hock and Schittkowski collection
% Initial guess: x0 = (2,2); f(x0) = log(5)-2 (not feasible)
% Optimal sol:   x* = (0,sqrt(3));  f(x*) = -sqrt(3)
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = (1+x(1)^2)^2 + x(2)^2 - 4;

end
