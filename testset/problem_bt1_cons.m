function c = problem_bt1_cons( x )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: problem 1 from P.T. Boggs and J.W. Tolle,
% "A strategy for global convergence in a sequential 
% quadratic programming algorithm", SINUM 26(3), pp. 600-623, 1989.
% Initial guess: x0 = (0.08,0.06)
% Optimal sol:   x* = (1,0)
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = x(1)^2+x(2)^2-1.0;

end
