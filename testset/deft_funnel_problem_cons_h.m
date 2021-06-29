function cons = deft_funnel_problem_cons_h(x, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Function called by 'run_deft_funnel.m' and used for running DEFT-FUNNEL
% on a collection of test problems. It sets the white-box constraint functions 
% of the test problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (size(varargin, 2) > 0)
   nprob = varargin{1};
else
   nprob = 1;
end

if (nprob == 1)
    
   % Problem: HS6
   cons = [];
   
elseif (nprob == 2)
    
   % Problem: HS7
   cons = [];
   
elseif (nprob == 3)
    
   % Problem: HS9
   cons = [];
   
elseif (nprob == 4)
    
   % Problem: HS10
   cons = [];
   
elseif (nprob == 5)
    
   % Problem: HS12
   cons = [];
   
elseif (nprob == 6)
    
   % Problem: HS13
   cons = [];
   
elseif (nprob == 7) 
   
   % Problem: HS21
   cons = [];

elseif (nprob == 8)
   
   % Problem: HS23
   cons = [];
   
elseif (nprob == 9)
   
   % Problem: BT1
   cons = [];
   
elseif (nprob == 10)
   
   % Problem: WB4
   cons = [];
   
elseif (nprob == 11)
   
   % Problem: GTCD4
   cons = [];
   
elseif (nprob == 12)
   
   % Problem: PVD4
   cons = [];
   
elseif (nprob == 13)
   
   % Problem: SR7
   cons = [];
   
elseif (nprob == 14)
   
   % Problem: Hesse
   cons = [];
   
elseif (nprob == 15)
   
   % Problem: Gomez #3
   cons = [];
   
elseif (nprob == 16)
   
   % Problem: G3
   cons = [];
   
elseif (nprob == 17)
   
   % Problem: G4
   cons = [];
   
elseif (nprob == 18)
   
   % Problem: G6
   cons = [];
   
elseif (nprob == 19)
   
   % Problem: G7
   cons = [];
   
elseif (nprob == 20)
   
   % Problem: G8
   cons = [];
   
elseif (nprob == 21)
   
   % Problem: G9
   cons = [];
   
elseif (nprob == 22)
   
   % Problem: G11
   cons = [];
   
elseif (nprob == 23)
   
   % Problem: Harley Pooling Problem
   cons = [];
   
elseif (nprob == 24)
   
   % Problem: GTCD4 (grey box)
   cons = problem_greybox_GTCD4_cons_h(x);
   
elseif (nprob == 25)
   
   % Problem: SR7 (grey box)
   cons = problem_greybox_SR7_cons_h(x);
   
elseif (nprob == 26)
   
   % Problem: Hesse (grey box)
   cons = problem_greybox_hesse_cons_h(x);
   
elseif (nprob == 27)
   
   % Problem: HS21 (grey box)
   cons = [];
   
elseif (nprob == 28)
   
   % Problem: HS23 (grey box)
   cons = problem_greybox_hs23_cons_h(x);
   
elseif (nprob == 29)
   
   % Problem: WB4 (grey box)
   cons = problem_greybox_WB4_cons_h(x);
   
end
