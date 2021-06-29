function f = problem_G16_obj(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: Problem G16 in  "Empirical analysis of a modified artificial bee 
% colony for constrained numerical optimization", Applied Mathematics and 
% Computation, 218(22):10943–10973, 2012,
% by E. Mezura-Montes and O. Cetina-Domíngue.
%
% Desc: 
%     - Number of variables: 5
%     - Number of constraints (not bounds): 38 inequalities
%     - Objective function: non-linear
%     - Constraints: non-linear
%
% Lower and upper bounds for the constraint(s):
% lc = -Inf, i=1,...,38
% uc = 0, i=1,...,38
%
% Lower and upper bounds for the decision variables x:
% lx = (704.4148, 68.6, 0, 193, 25)
% ux = (906.3855, 288.88, 134.75, 287.0966, 84.1988)
%
% Global optimal sol (best known): 
% x* = (705.174537070090537, 68.5999999999999943, 102.899999999999991, ...
%       282.324931593660324, 37.5841164258054832)
% f(x*) = -1.90515525853479
%
% Programming: Phillipe R. Sampaio
% This file is part of the DEFT-FUNNEL software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y1 = x(2) + x(3) + 41.6;
c1 = 0.024*x(4) - 4.62;
y2 = 12.5/c1 + 12;
c2 = 0.0003535*x(1)^2 + 0.5311*x(1) + 0.08705*y2*x(1);
c3 = 0.052*x(1) + 78 + 0.002377*y2*x(1);
y3 = c2/c3;
y4 = 19*y3;
c4 = 0.04782*(x(1)-y3) + (0.1956*(x(1)-y3)^2)/x(2) + 0.6376*y4 + 1.594*y3;
c5 = 100*x(2);
c6 = x(1) - y3 - y4;
c7 = 0.950 - c4/c5;
y5 = c6*c7;
y6 = x(1) - y5 - y4 - y3;
c8 = (y5 + y4)*0.995;
y7 = c8/y1;
y8 = c8/3798;
c9 = y7 - 0.0663*y7/y8 - 0.3153;
y9 = 96.82/c9 + 0.321*y1;
y10 = 1.29*y5 + 1.258*y4 + 2.29*y3 + 1.71*y6;
y11 = 1.71*x(1) - 0.452*y4 + 0.580*y3;
c10 = 12.3/752.3;
c11 = (1.75*y2)*(0.995*x(1));
c12 = 0.995*y10 + 1998;
y12 = c10*x(1) + c11/c12;
y13 = c12 + 1.75*y2;
y14 = 3623 + 64.4*x(2) + 58.4*x(3) + 146312/(y9 + x(5));
c13 = 0.995*y10 + 60.8*x(2) + 48*x(4) - 0.1121*y14 - 5095;
y15 = y13/c13;
y16 = 148000 - 331000*y15 + 40*y13 - 61*y15*y13;
c14 = 2324*y10 - 28740000*y2;
y17 = 14130000 - 1328*y10 - 531*y11 + c14/c12;
c15 = y13/y15 - y13/0.52;
c16 = 1.104 - 0.72*y15;
c17 = y9 + x(5);

f = 0.000117*y14 + 0.1365 + 0.00002358*y13 + 0.000001502*y16 + 0.0321*y12 + ...
    0.004324*y5 + 0.0001*(c15/c16) + 37.48*y2/c12 - 0.0000005843*y17;

end