function c = problem_hs9_cons( x )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem 9 from Hock and Schittkowski collection
% Initial guess: x0 = (0,0); f(x0) = 0 (feasible)
% Optimal sol:   x* = (12k-3, 16k-4), k=0, +-1, +-2,... ;  f(x*) = -0.5
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 4*x(1) - 3*x(2);

end
