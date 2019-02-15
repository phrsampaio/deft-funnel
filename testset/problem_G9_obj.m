function f = problem_G9_obj(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem G9 in  "Evolutionary algorithms for constrained parameter 
% optimization problems.", Evolutionary Computation, 4(1):1–32, 1996,
% by Michalewicz, Z., Schoenauer, M.
%
% Desc: 
%     - Number of variables: 7
%     - Number of constraints (not bounds): 4
%     - Objective function: non-linear
%     - Constraints: non-linear
%     - Type: First and last constraints active in the global optimum
%
% Global optimal sol:   
% x*    = (2.330499, 1.951372, -0.4775414, 4.365726, -0.6244870, 1.038131, 1.594227)
% f(x*) = 680.6300573
% lb    = -10, i=1..7
% ub    = 10, i=1..7
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = (x(1)-10)^2 + 5*(x(2)-12)^2 + x(3)^4 + 3*(x(4)-11)^2 + ...
    10*x(5)^6 + 7*x(6)^2 + x(7)^4 - 4*x(6)*x(7) - 10*x(6) - 8*x(7);

end