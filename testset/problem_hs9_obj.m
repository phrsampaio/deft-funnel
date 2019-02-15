function f = problem_hs9_obj( x ) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem 9 from Hock and Schittkowski collection
% Initial guess: x0 = (0,0); f(x0) = 0 (feasible)
% Optimal sol:   x* = (12k-3, 16k-4), k=0, +-1, +-2,... ;  f(x*) = -0.5
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = sin(pi*x(1)/12)*cos(pi*x(2)/16);

end