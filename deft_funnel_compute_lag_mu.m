function [ mu, ws, zs, wx, zx ] = deft_funnel_compute_lag_mu( iterate, derivatives, setting )
                                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Computes estimates for the Lagrange multipliers using the solver 
% BLLS (bound-constrained linear least squares) implemented in
% 'deft_funnel_blls_exp.m' and 'deft_funnel_blls_spwmin.m'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Retrieve necessary data
g = derivatives.gfx;
J = derivatives.J;
x = iterate.x;
s = iterate.s;

ls = setting.ls;
us = setting.us;
lx = setting.lx( iterate.indfree );
ux = setting.ux( iterate.indfree );

n = length( x );
m = size( J, 1 );
Im = eye( m );
In = eye( n );
size_mu = m;

% Define tolerance for active bounds
tol = 1.0e-12;

% Consider only the active upper s bounds in the computation
indws = find( abs(s-us) < tol );
Iws = Im( :, indws );
size_wsAct = size( Iws, 2 );

% Consider only the active lower s bounds in the computation
indzs = find( abs(ls-s) < tol );
Izs = -Im( :, indzs );
size_zsAct = size( Izs, 2 );

% Consider only the active upper x bounds in the computation
indwx = find( abs(x-ux) < tol );
Iwx = In( :, indwx );
size_wxAct = size( Iwx, 2 );

% Consider only the active lower x bounds in the computation
indzx = find( abs(lx-x) < tol );
Izx = -In( :, indzx );
size_zxAct = size( Izx, 2 );

zerosn2m = zeros( n, size( Iws, 2 ) + size( Izs, 2 ) );
zerosm2n = zeros( m, size( Iwx, 2 ) + size( Izx, 2 ) );

glagMatrix = [ J' zerosn2m Iwx Izx; -Im Iws Izs zerosm2n ];
size_glagMatrix = size( glagMatrix, 2 );

b = [ g; zeros( m, 1 ) ];

lb( 1:size_mu ) = -Inf;
ub( 1:size_mu ) = Inf;
lb( size_mu+1:size_glagMatrix ) = 0;
ub( size_mu+1:size_glagMatrix ) = Inf;

% Solve the bound-constrained least-squares problem to obtain the 
% Lagrangian multipliers
[ lagmu_temp, resn, exitc ] = deft_funnel_blls_exp( glagMatrix, ...
    -b, lb',ub' );

mu = lagmu_temp( 1:size_mu );

ws = zeros( m, 1 );
if size_wsAct > 0
    ws( indws ) = lagmu_temp( size_mu+1:size_mu+size_wsAct );
end

zs = zeros(m,1);
if size_zsAct > 0
    zs( indzs ) = lagmu_temp( size_mu+size_wsAct+1:size_mu+size_wsAct+size_zsAct );
end

wx = zeros(n,1);
if size_wxAct > 0
    wx( indwx ) = lagmu_temp( size_mu+size_wsAct+size_zsAct+1:size_mu+size_wsAct+size_zsAct+size_wxAct );
end

zx = zeros(n,1);
if size_zxAct > 0
    zx( indzx ) = lagmu_temp( size_mu+size_wsAct+size_zsAct+size_wxAct+1:size_glagMatrix );
end

end % end of deft_funnel_compute_lag_mu
