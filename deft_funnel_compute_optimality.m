function [ pi_f, exitc ] = deft_funnel_compute_optimality( iterate, derivatives, setting )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Computes the dual optimiality measure by solving the following problem
%
%       min          g_n * dr
%       s.t.      J_s * dr = 0,
%             lx - x <= dr^x <= ux - x,
%             ls - s <= dr^s <= us - s,
%                 ||dr|| < = 1.
%
% The optimality measure if defined as
%              pi_f = | <g_n,dr*> |, 
% where dr* is the solution to the problem above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exitc = 0; % Initialize exit condition to successful
n = length( iterate.x );
m = length( iterate.s );
J_s = derivatives.J_s;
g_n = [ derivatives.gfx; zeros( m, 1 ) ];

lb(1:n+m) = -1.0;
ub(1:n+m) = 1.0;

dlx = setting.lx( iterate.indfree ) - iterate.x;
dlx = dlx';
dux = setting.ux( iterate.indfree ) - iterate.x;
dux = dux';
lb(1:n) = max( lb(1:n), dlx );
ub(1:n) = min( ub(1:n), dux );

dls = setting.ls - iterate.s;
dls = dls';
dus = setting.us - iterate.s;
dus = dus';
lb(n+1:n+m) = min( max( lb( n+1:n+m ), dls ), setting.us' );
ub(n+1:n+m) = max( min( ub( n+1:n+m ), dus ), setting.ls' );

initPoint = zeros( n+m, 1 );
b = zeros( m, 1 );

options = optimset( 'LargeScale','off','Simplex','off','Display','off' );
[ dr, fevallinp, exitflag, output ] = linprog( g_n, [], [], J_s, b, ...
    lb', ub', initPoint, options );

% Compute the dual optimality measure
if ( norm( dr ) > 1.0e-14 || exitflag == -2 )
    
    pi_f = abs(g_n.' * dr);

    if ( exitflag == -2 ) % Try another method

        options = optimset( 'LargeScale','off','Simplex','on',     ...
            'Display','off' );
        [ dr, fevallinp, exitflag, output ] = linprog(g_n, [], [],  ...
            J_s, b, lb', ub', initPoint, options );
        pi_f = abs(g_n.' * dr);

        if ( exitflag == -2 ) % bad news
            exitc = 1;
        end
    end
else
    pi_f = 0.0;
end

end % end of deft_funnel_compute_optimality