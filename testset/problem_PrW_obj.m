function f = problem_PrW_obj(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Welded Beam Design Problem in "An efficient constraint handling 
% method for genetic algorithms", Computer Methods in Applied Mechanics and 
% Engineering 186, 311–338, 2002, by Deb, K.
%
% Desc: 
%     - Number of variables: 4
%     - Number of constraints (not bounds): 4
%     - Objective function: non-linear
%     - Constraints: non-linear
%     - Type: engineering design problem
%
% Global optimal sol obtained by the FSA algorithm: 
% x*    = (0.20564426101885, 3.47257874213172, 9.03662391018928, 0.20572963979791)
% f(x*) = 1.728226 (new:  1.7250022 by FSA)
% lb    = (0.125, 0.1, 0.1, 0.1)
% ub    = (10, 10, 10, 10)
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = 1.10471*x(1)^2*x(2) + 0.04811*x(3)*x(4)*(14.0 + x(2));

end