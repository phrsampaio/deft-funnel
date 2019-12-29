function [f, g] = deft_funnel_tangent_prob(t, M, g_n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Evaluates the objective function and the gradient of the tangent
% step problem at the point 't' given as input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = g_n.'*t + 0.5 * t.'* (M * t);
g = g_n + M * t;
  
end
