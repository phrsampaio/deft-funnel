function cons = deft_funnel_problem_cons_c(x, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Function called by 'run_deft_funnel.m' and used for running DEFT-FUNNEL
% on a collection of test problems. It sets the black-box constraint functions 
% of the test problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (size(varargin, 2) > 0)
   nprob = varargin{1};
else
   nprob = 1;
end

if (nprob == 1)
    
   % Problem: HS6
   cons = problem_hs6_cons(x);
   
elseif (nprob == 2)
    
   % Problem: HS7
   cons = problem_hs7_cons(x);
   
elseif (nprob == 3)
    
   % Problem: HS9
   cons = problem_hs9_cons(x);
   
elseif (nprob == 4)
    
   % Problem: HS10
   cons = problem_hs10_cons(x);
   
elseif (nprob == 5)
    
   % Problem: HS12
   cons = problem_hs12_cons(x);
   
elseif (nprob == 6)
    
   % Problem: HS13
   cons = problem_hs13_cons(x);
   
elseif (nprob == 7) 
   
   % Problem: HS21
   cons = problem_hs21_cons(x);

elseif (nprob == 8)
   
   % Problem: HS23
   cons = problem_hs23_cons(x);
   
elseif (nprob == 9)
   
   % Problem: BT1
   cons = problem_bt1_cons(x);
   
elseif (nprob == 10)
   
   % Problem: WB4
   cons = problem_WB4_cons(x);
   
elseif (nprob == 11)
   
   % Problem: GTCD4
   cons = problem_GTCD4_cons(x);
   
elseif (nprob == 12)
   
   % Problem: PVD4
   cons = problem_PVD4_cons(x);
   
elseif (nprob == 13)
   
   % Problem: SR7
   cons = problem_SR7_cons(x);
   
elseif (nprob == 14)
   
   % Problem: Hesse
   cons = problem_hesse_cons(x);
   
elseif (nprob == 15)
   
   % Problem: Gomez #3
   cons = problem_gomez_pb3_cons(x);
   
elseif (nprob == 16)
   
   % Problem: G3
   cons = problem_G3_cons(x);
   
elseif (nprob == 17)
   
   % Problem: G4
   cons = problem_G4_cons(x);
   
elseif (nprob == 18)
   
   % Problem: G6
   cons = problem_G6_cons(x);
   
elseif (nprob == 19)
   
   % Problem: G7
   cons = problem_G7_cons(x);
   
elseif (nprob == 20)
   
   % Problem: G8
   cons = problem_G8_cons(x);
   
elseif (nprob == 21)
   
   % Problem: G9
   cons = problem_G9_cons(x);
   
elseif (nprob == 22)
   
   % Problem: G11
   cons = problem_G11_cons(x);
   
elseif (nprob == 23)
   
   % Problem: Harley Pooling Problem
   cons = problem_harley_cons(x);
   
elseif (nprob == 24)
   
   % Problem: GTCD4 (grey box)
   cons = [];
   
elseif (nprob == 25)
   
   % Problem: SR7 (grey box)
   cons = problem_greybox_SR7_cons_c(x);
   
elseif (nprob == 26)
   
   % Problem: Hesse (grey box)
   cons = problem_greybox_hesse_cons_c(x);
   
elseif (nprob == 27)
   
   % Problem: HS21 (grey box)
   cons = problem_greybox_hs21_cons_c(x);
   
elseif (nprob == 28)
   
   % Problem: HS23 (grey box)
   cons = problem_greybox_hs23_cons_c(x);
   
end
