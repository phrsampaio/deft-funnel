function c = problem_hs23_cons( x )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem 23 from Hock and Schittkowski collection
% Initial guess: x0 = (3,1); f(x0) = 10 (not feasible)
% Optimal sol:   x* = (1,1);  f(x*) = 2
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c(1) = x(1) + x(2) - 1;
c(2) = x(1)^2 + x(2)^2 - 1;
c(3) = 9*x(1)^2 + x(2)^2 - 9;
c(4) = x(1)^2 - x(2);
c(5) = x(2)^2 - x(1);
   
end
