function h = problem_greybox_GTCD4_cons_h(x)

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

% White-box constraint
h = x(4)/x(2) + 1/x(2) - 1;

end
