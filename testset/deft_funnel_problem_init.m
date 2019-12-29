function [x0, n, nb_cons_c, nb_cons_h, ls, us, lx, ux, type_f] = deft_funnel_problem_init(nprob)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Function called by 'run_deft_funnel.m' and used for running DEFT-FUNNEL
% on a collection of test problems. It sets the entry parameters of the
% test problem and defines if a constraint in 'deft_funnel_problem_cons.m'
% is an equality or an inequality through the lower bounds 'ls' and 
% the upper bounds 'us'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nprob == 1)

    % Problem: HS6
    x0        = [-1.2 1];
    n         = 2;
    nb_cons_c = 1;
    nb_cons_h = 0;
    ls        = 0;
    us        = 0;
    lx        = [-Inf -Inf];
    ux        = [Inf Inf];
    type_f    = 'BB';

elseif (nprob == 2)

    % Problem: HS7
    x0        = [2 2];
    n         = 2;
    nb_cons_c = 1;
    nb_cons_h = 0;
    ls        = 0;
    us        = 0;
    lx        = [-Inf -Inf];
    ux        = [Inf Inf];
    type_f    = 'BB';

elseif (nprob == 3)

    % Problem: HS9
    x0        = [0 0];
    n         = 2;
    nb_cons_c = 1;
    nb_cons_h = 0;
    ls        = 0;
    us        = 0;
    lx        = [-Inf -Inf];
    ux        = [Inf Inf];
    type_f    = 'BB';

elseif (nprob == 4)

    % Problem: HS10
    x0        = [-10 10];
    n         = 2;
    nb_cons_c = 1;
    nb_cons_h = 0;
    ls        = 0;
    us        = Inf;
    lx        = [-Inf -Inf];
    ux        = [Inf Inf];
    type_f    = 'BB';

elseif (nprob == 5)

    % Problem: HS12
    x0        = [0 0];
    n         = 2;
    nb_cons_c = 1;
    nb_cons_h = 0;
    ls        = 0;
    us        = Inf;
    lx        = [-Inf -Inf];
    ux        = [Inf Inf];
    type_f    = 'BB';

elseif (nprob == 6)

    % Problem: HS13
    x0        = [-2 -2];
    n         = 2;
    nb_cons_c = 1;
    nb_cons_h = 0;
    ls        = 0;
    us        = Inf;
    lx        = [0 0];
    ux        = [Inf Inf];
    type_f    = 'BB';

elseif (nprob == 7)

    % Problem: HS21
    x0        = [-1 -1];
    n         = 2;
    nb_cons_c = 1;
    nb_cons_h = 0;
    ls        = 0;
    us        = Inf;
    lx        = [2 -50];
    ux        = [50 50];
    type_f    = 'BB';

elseif (nprob == 8)

    % Problem: HS23
    x0        = [3 1];
    n         = 2;
    nb_cons_c = 5;
    nb_cons_h = 0;
    ls        = [0 0 0 0 0];
    us        = [Inf Inf Inf Inf Inf];
    lx        = [-50 -50];
    ux        = [50 50];
    type_f    = 'BB';

elseif (nprob == 9)

    % Problem: BT1
    x0        = [0.08 0.06];
    n         = 2;
    nb_cons_c = 1;
    nb_cons_h = 0;
    ls        = 0;
    us        = 0;
    lx        = [-Inf -Inf];
    ux        = [Inf Inf];
    type_f    = 'BB';

elseif (nprob == 10)

    % Problem: WB4
    x0        = [1 1 1 1];
    n         = 4;
    nb_cons_c = 6;
    nb_cons_h = 0;
    ls        = [-Inf -Inf -Inf -Inf -Inf -Inf];
    us        = [0 0 0 0 0 0];
    lx        = [0.125 0.1 0.1 0.1];
    ux        = [10 10 10 10];
    type_f    = 'BB';

elseif (nprob == 11)

    % Problem: GTCD4
    x0        = [20 5 30 20];
    n         = 4;
    nb_cons_c = 1;
    nb_cons_h = 0;
    ls        = -Inf;
    us        = 0;
    lx        = [20 1 20 0.1];
    ux        = [50 10 50 60];
    type_f    = 'BB';

elseif (nprob == 12)

    % Problem: PVD4
    x0        = [0 0 10 100];
    n         = 4;
    nb_cons_c = 3;
    nb_cons_h = 0;
    ls        = [-Inf -Inf -Inf];
    us        = [0 0 0];
    lx        = [0 0 0 0];
    ux        = [1 1 50 240];
    type_f    = 'BB';

elseif (nprob == 13)

    % Problem: SR7
    x0        = [2.6 0.7 17 7.3 7.3 2.9 5];
    n         = 7;
    nb_cons_c = 11;
    nb_cons_h = 0;
    ls        = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf];
    us        = [0 0 0 0 0 0 0 0 0 0 0];
    lx        = [2.6 0.7 17 7.3 7.3 2.9 5];
    ux        = [3.6 0.8 28 8.3 8.3 3.9 5.5];
    type_f    = 'BB';

elseif (nprob == 14)

    % Problem: Hesse
    %x0        = [1 1 1 0 1 0];
    x0        = [0.177468345583169 0.400230345036942 4.93626491163683 0.0650223791902058 3.82581226100484 7.18570456709047];
    n         = 6;
    nb_cons_c = 5;
    nb_cons_h = 0;
    ls        = [4 4 -Inf -Inf 2];
    us        = [Inf Inf 2 2 6];
    lx        = [0 0 1 0 1 0];
    ux        = [Inf Inf 5 6 5 10];
    type_f    = 'BB';

elseif (nprob == 15)

    % Problem: Gomez #3
    x0        = [-1 0];
    n         = 2;
    nb_cons_c = 1;
    nb_cons_h = 0;
    ls        = -Inf;
    us        = 0;
    lx        = [-1 -1];
    ux        = [1 1];
    type_f    = 'BB';

elseif (nprob == 16)

    % Problem: G3
    x0        = [0 0];
    n         = 2;
    nb_cons_c = 1;
    nb_cons_h = 0;
    ls        = 0;
    us        = 0;
    lx        = [0 0];
    ux        = [1 1];
    type_f    = 'BB';
   
elseif (nprob == 17)

    % Problem: G4
    x0        = [78 33 27 27 27];
    n         = 5;
    nb_cons_c = 3;
    nb_cons_h = 0;
    ls        = [0 90 20];
    us        = [92 110 25];
    lx        = [78 33 27 27 27];
    ux        = [102 45 45 45 45];
    type_f    = 'BB';
   
elseif (nprob == 18)

    % Problem: G6
    x0        = [20.1, 5.84];
    n         = 2;
    nb_cons_c = 2;
    nb_cons_h = 0;
    ls        = [0 0];
    us        = [Inf Inf];
    lx        = [13 0];
    ux        = [100 100];
    type_f    = 'BB';
   
elseif (nprob == 19)

    % Problem: G7
    x0        = [0 0 0 0 0 0 0 0 0 0];
    n         = 10;
    nb_cons_c = 8;
    nb_cons_h = 0;
    ls        = [0 0 0 0 0 0 0 0];
    us        = [Inf Inf Inf Inf Inf Inf Inf Inf];
    lx        = [-10 -10 -10 -10 -10 -10 -10 -10 -10 -10];
    ux        = [10 10 10 10 10 10 10 10 10 10];
    type_f    = 'BB';
   
elseif (nprob == 20)

    % Problem: G8
    x0        = [0 0];
    n         = 2;
    nb_cons_c = 2;
    nb_cons_h = 0;
    ls        = [-Inf -Inf];
    us        = [0 0];
    lx        = [0 0];
    ux        = [10 10];
    type_f    = 'BB';
   
elseif (nprob == 21)

    % Problem: G9
    x0        = [0 0 0 0 0 0 0];
    n         = 7;
    nb_cons_c = 4;
    nb_cons_h = 0;
    ls        = [0 0 0 0];
    us        = [Inf Inf Inf Inf];
    lx        = [-10 -10 -10 -10 -10 -10 -10];
    ux        = [10 10 10 10 10 10 10];
    type_f    = 'BB';
   
elseif (nprob == 22)

    % Problem: G11
    x0        = [0 0];
    n         = 2;
    nb_cons_c = 1;
    nb_cons_h = 0;
    ls        = 0;
    us        = 0;
    lx        = [-1 -1];
    ux        = [1 1];
    type_f    = 'BB';
   
elseif (nprob == 23)

    % Problem: Harley pooling problem
    x0        = [0 0 0 0 0 0 0 0 0];
    n         = 9;
    nb_cons_c = 6;
    nb_cons_h = 0;
    ls        = [0 0 0 -Inf -Inf 0];
    us        = [0 0 0 0 0 0];
    lx        = [0 0 0 0 0 0 0 0 0];
    ux        = [600 200 500 500 500 500 500 500 500];
    type_f    = 'BB';
   
elseif (nprob == 24)

    % Problem: GTCD4
    x0        = [20 5 30 20];
    n         = 4;
    nb_cons_c = 0;
    nb_cons_h = 1;
    ls        = -Inf;
    us        = 0;
    lx        = [20 1 20 0.1];
    ux        = [50 10 50 60];
    type_f    = 'BB';
   
elseif (nprob == 25)

    % Problem: SR7
    x0        = [2.6 0.7 17 7.3 7.3 2.9 5];
    n         = 7;
    nb_cons_c = 9;
    nb_cons_h = 2;
    ls        = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf];
    us        = [0 0 0 0 0 0 0 0 0 0 0];
    lx        = [2.6 0.7 17 7.3 7.3 2.9 5];
    ux        = [3.6 0.8 28 8.3 8.3 3.9 5.5];
    type_f    = 'WB';
   
elseif (nprob == 26)

    % Problem: Hesse
    %x0        = [1 1 1 0 1 0];
    x0        = [0.177468345583169 0.400230345036942 4.93626491163683 0.0650223791902058 3.82581226100484 7.18570456709047];
    n         = 6;
    nb_cons_c = 3;
    nb_cons_h = 2;
    ls        = [-Inf -Inf 2 4 4];
    us        = [2 2 6 Inf Inf];
    lx        = [0 0 1 0 1 0];
    ux        = [Inf Inf 5 6 5 10];
    type_f    = 'WB';
   
elseif (nprob == 27)

    % Problem: HS21
    x0        = [-1 -1];
    n         = 2;
    nb_cons_c = 1;
    nb_cons_h = 0;
    ls        = 0.0;
    us        = Inf;
    lx        = [2 -50];
    ux        = [50 50];
    type_f    = 'WB';
   
elseif (nprob == 28)

    % Problem: HS23
    x0        = [3 1];
    n         = 2;
    nb_cons_c = 2;
    nb_cons_h = 3;
    ls        = [0 0 0 0 0];
    us        = [Inf Inf Inf Inf Inf];
    lx        = [-50 -50];
    ux        = [50 50];
    type_f    = 'WB';
   
end % end of deft_funnel_init
