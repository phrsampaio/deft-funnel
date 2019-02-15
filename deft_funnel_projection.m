function y = deft_funnel_projection( x, Aeq, beq, lb, ub )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Solves the quadratic convex problem below:
%
%    min_y     ||y - x||^2
%    s.t.      Aeq*y = Beq, 
%            lb <= y <= ub.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = eye( length( x ) );
f = -x;

options.LargeScale = 'off';
options.Display = 'none';

[ y, feval, exitflag ] = quadprog( H, f, [], [], ...
    Aeq, beq, lb, ub, x, options);

if ( exitflag == -2 )
    disp(' ');
    disp(' Projection problem infeasible');
    disp(' Setting tangent step to zero');
    disp(' ');
    y = zeros( length( x ), 1 );
end

end % end of deft_funnel_projection