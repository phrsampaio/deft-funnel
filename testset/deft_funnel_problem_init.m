function [ x0, nbcons, ls, us, lx, ux ] = deft_funnel_problem_init( nprob )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Function called by 'run_deft_funnel.m' and used for running DEFT-FUNNEL
% on a collection of test problems. It sets the entry parameters of the
% test problem and defines if a constraint in 'deft_funnel_problem_cons.m'
% is an equality of an inequality through the lower bounds 'ls' and 
% the upper bounds 'us'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( nprob == 1 )

    % Problem: HS6
    x0      = [ -1.2 1 ];
    nbcons  = 1;
    ls      = 0;
    us      = 0;
    lx      = [ -Inf -Inf ];
    ux      = [ Inf Inf ];

elseif ( nprob == 2 )

    % Problem: HS7
    x0      = [ 2 2 ];
    nbcons  = 1;
    ls      = 0;
    us      = 0;
    lx      = [ -Inf -Inf ];
    ux      = [ Inf Inf ];

elseif ( nprob == 3 )

    % Problem: HS9
    x0     = [ 0 0 ];
    nbcons = 1;
    ls     = 0;
    us     = 0;
    lx     = [ -Inf -Inf ];
    ux     = [ Inf Inf ];

elseif ( nprob == 4 )

    % Problem: HS10
    x0     = [ -10 10 ];
    nbcons = 1;
    ls     = 0;
    us     = Inf;
    lx     = [ -Inf -Inf ];
    ux     = [ Inf Inf ];

elseif ( nprob == 5 )

    % Problem: HS12
    x0     = [ 0 0 ];
    nbcons = 1;
    ls     = 0;
    us     = Inf;
    lx     = [ -Inf -Inf ];
    ux     = [ Inf Inf ];

elseif ( nprob == 6 )

    % Problem: HS13
    x0     = [ -2 -2 ];
    nbcons = 1;
    ls     = 0;
    us     = Inf;
    lx     = [ 0 0 ];
    ux     = [ Inf Inf ];

elseif ( nprob == 7 )

    % Problem: HS21
    x0     = [ -1 -1 ];
    nbcons = 1;
    ls     = 0.0;
    us     = Inf;
    lx     = [ 2 -50 ];
    ux     = [ 50 50 ];

elseif ( nprob == 8 )

    % Problem: HS23
    x0     = [ 3 1 ];
    nbcons = 5;
    ls     = [ 0 0 0 0 0 ];
    us     = [ Inf Inf Inf Inf Inf ];
    lx     = [ -50 -50 ];
    ux     = [ 50 50 ];

elseif ( nprob == 9 )

    % Problem: BT1
    x0     = [0.08 0.06];
    nbcons = 1;
    ls     = 0;
    us     = 0;
    lx     = [ -Inf -Inf ];
    ux     = [ Inf Inf ];

elseif ( nprob == 10 )

    % Problem: Test problem 3 in "Quadratically Constrained Problems" section of
    % "Handbook of test problems in local and Global Optimization"
    %x0     = [1 1 1 0 1 0];
    x0     = [ 0.177468345583169 0.400230345036942 4.93626491163683 0.0650223791902058 3.82581226100484 7.18570456709047 ];
    nbcons = 5;
    ls     = [ 4 4 -Inf -Inf 2];
    us     = [ Inf Inf 2 2 6 ];
    lx     = [ 0 0 1 0 1 0 ];
    ux     = [ Inf Inf 5 6 5 10 ];

end
   
end % end of deft_funnel_init
