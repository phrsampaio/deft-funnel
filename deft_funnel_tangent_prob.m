function [ f, g ] = deft_funnel_tangent_prob( t, M, g_n )

f = g_n.'*t + 0.5 * t.'* ( M * t );
g = g_n + M * t;
  
end
