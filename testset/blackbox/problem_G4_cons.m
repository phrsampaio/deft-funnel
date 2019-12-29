function c = problem_G4_cons(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem G4 in  "Evolutionary algorithms for constrained parameter 
% optimization problems.", Evolutionary Computation, 4(1):1–32, 1996,
% by Michalewicz, Z., Schoenauer, M.
%
% Desc: 
%     - Number of variables: 5
%     - Number of constraints (not bounds): 6 (3 double inequalities)
%     - Objective function: non-linear
%     - Constraints: non-linear
%     - Two constraints (upper bound of the frst inequality and the lower
%       bound of the third inequality) are active at the optimum.
%
% Lower and upper bounds for the constraint(s):
% lc = (0, 90, 20)
% uc = (92, 110, 25)
%
% % Lower and upper bounds for the decision variables x:
% lx = (78, 33, 27, 27, 27)
% ux = (102, 45, 45, 45, 45)
%
% Global optimal sol: x* = (78, 33, 29.995, 45, 36.776); f(x*) = -30665.539
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c(1) = 85.334407 + 0.0056858*x(2)*x(5) + 0.00026*x(1)*x(4) - 0.0022053*x(3)*x(5);
c(2) = 80.51249 + 0.0071317*x(2)*x(5) + 0.0029955*x(1)*x(2) + 0.0021813*x(3)^2;
c(3) = 9.300961 + 0.0047026*x(3)*x(5) + 0.0012547*x(1)*x(3) + 0.0019085*x(3)*x(4);

end