function f = deft_funnel_problem_obj(x, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Function called by 'run_deft_funnel.m' and used for running DEFT-FUNNEL
% on a collection of test problems. It sets the objective function of the
% test problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (size(varargin, 2) > 0)
   nprob = varargin{1};
else
   nprob = 1;
end

if (nprob == 1)
    
   % Problem: HS6
   f = problem_hs6_obj(x);
   
elseif (nprob == 2)
    
   % Problem: HS7
   f = problem_hs7_obj(x);
   
elseif (nprob == 3)
    
   % Problem: HS9
   f = problem_hs9_obj(x);
   
elseif (nprob == 4)
    
   % Problem: HS10
   f = problem_hs10_obj(x);
   
elseif (nprob == 5)
    
   % Problem: HS12
   f = problem_hs12_obj(x);
   
elseif (nprob == 6)
    
   % Problem: HS13
   f = problem_hs13_obj(x);
   
elseif (nprob == 7) 
   
   % Problem: HS21
   f = problem_hs21_obj(x);

elseif (nprob == 8)
    
   % Problem: HS23
   f = problem_hs23_obj(x);
   
elseif (nprob == 9)
    
   % Problem: BT1
   f = problem_bt1_obj(x);
   
elseif (nprob == 10)
   
   % Problem: WB4
   f = problem_WB4_obj(x);
   
elseif (nprob == 11)
   
   % Problem: GTCD4
   f = problem_GTCD4_obj(x);
   
elseif (nprob == 12)
   
   % Problem: PVD4
   f = problem_PVD4_obj(x);
   
elseif (nprob == 13)
   
   % Problem: SR7
   f = problem_SR7_obj(x);
   
elseif (nprob == 14)
   
   % Problem: Hesse
   f = problem_hesse_obj(x);
   
elseif (nprob == 15)
   
   % Problem: Gomez #3
   f = problem_gomez_pb3_obj(x);
   
elseif (nprob == 16)
   
   % Problem: G3
   f = problem_G3_obj(x);
   
elseif (nprob == 17)
   
   % Problem: G4
   f = problem_G4_obj(x);
   
elseif (nprob == 18)
   
   % Problem: G6
   f = problem_G6_obj(x);
   
elseif (nprob == 19)
   
   % Problem: G7
   f = problem_G7_obj(x);
   
elseif (nprob == 20)
   
   % Problem: G8
   f = problem_G8_obj(x);
   
elseif (nprob == 21)
   
   % Problem: G9
   f = problem_G9_obj(x);
   
elseif (nprob == 22)
   
   % Problem: G11
   f = problem_G11_obj(x);
   
elseif (nprob == 23)
   
   % Problem: Harley Pooling Problem
   f = problem_harley_obj(x);
   
elseif (nprob == 24)
   
   % Problem: GTCD4 (grey box)
   f = problem_greybox_GTCD4_obj(x);
   
elseif (nprob == 25)
   
   % Problem: SR7 (grey box)
   f = problem_greybox_SR7_obj(x);
   
elseif (nprob == 26)
   
   % Problem: Hesse (grey box)
   f = problem_greybox_hesse_obj(x);
   
elseif (nprob == 27)
   
   % Problem: HS21 (grey box)
   f = problem_greybox_hs21_obj(x);
   
elseif (nprob == 28)
   
   % Problem: HS23 (grey box)
   f = problem_greybox_hs23_obj(x);
   
end

