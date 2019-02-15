function t = deft_funnel_quad_interp( a, f0, fa, g0 )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Desc: Called by the function deft_funnel_spg.m, this functions computes a
% quadratic interpolation for finding a step size.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = ( g0 * a^2 )/( 2*( g0*a + f0 - fa ) );

end % end of deft_funnel_quad_interp