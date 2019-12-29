function f = problem_greybox_GTCD4_obj(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: "Applied Geometric Programming", John Wiley & Sons Inc, New York, 
% 1976, by C. S. Beightler and D. T. Phillips.
%
% Desc: 
%     - Number of variables: 4
%     - Number of constraints (not bounds): 1 inequality (white box)
%     - Objective function: non-linear (black-box)
%     - Constraints: non-linear
%
% Lower and upper bounds for the constraint(s):
% lc = -Inf
% uc = 0
%
% Lower and upper bounds for the decision variables x:
% lx = (20, 1, 20, 0.1)
% ux = (50, 10, 50, 60)
%
% Global optimal sol (best known): 
% f(x*) = 2964893.85
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = ((x(1)^(0.5))*x(2)) / ((x(3)^(2/3))*(x(4)^(1/2)));
B = (x(2)^(0.219))/x(1);

f = (8.61*100000)*A + (3.69*10000)*x(3) + (7.72*100000000)*B - (765.43*1000000)/x(1);

end
