function cons = deft_funnel_problem_cons( x, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Function called by 'run_deft_funnel.m' and used for running DEFT-FUNNEL
% on a collection of test problems. It sets the constraint functions of the
% test problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( size( varargin, 2 ) > 0 )
   nprob = varargin{ 1 };
else
   nprob = 1;
end

if ( nprob == 1 )
    
   % Problem: HS6
   cons = problem_hs6_cons( x );
   
elseif ( nprob == 2 )
    
   % Problem: HS7
   cons = problem_hs7_cons( x );
   
elseif ( nprob == 3 )
    
   % Problem: HS9
   cons = problem_hs9_cons( x );
   
elseif ( nprob == 4 )
    
   % Problem: HS10
   cons = problem_hs10_cons( x );
   
elseif ( nprob == 5 )
    
   % Problem: HS12
   cons = problem_hs12_cons( x );
   
elseif ( nprob == 6 )
    
   % Problem: HS13
   cons = problem_hs13_cons( x );
   
elseif ( nprob == 7 ) 
   
   % Problem: HS21
   cons = problem_hs21_cons( x );

elseif ( nprob == 8 )
   
   % Problem: HS23
   cons = problem_hs23_cons( x );
   
elseif ( nprob == 9 )
   
   % Problem: BT1
   cons = problem_bt1_cons( x );
   
elseif ( nprob == 10 )
   
   % Problem: Test problem 3 in "Quadratically Constrained Problems" section of
   % "Handbook of test problems in local and Global Optimization"
   cons = problem_handbook_quadcons_pb3_cons( x );
   
end
