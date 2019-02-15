function t = deft_funnel_cub_interp( a, b, fa, fb, ga, gb )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Desc: Called by the function deft_funnel_spg.m, this functions computes a 
% cubic interpolation for finding a step size.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = ga + gb + 3*( fa-fb )/( b-a );
w = sqrt( z^2 - ga*gb );
if isreal( w )
    temp = b - ( b - a )*( gb + w - z )/( gb - ga + 2*w );
    t = min( max( temp, a ), b );
else
    t = mean( [ a b ] );
end

end % end of deft_funnel_cub_interp