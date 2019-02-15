function f = deft_funnel_problem_obj( x, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function called by 'run_deft_funnel.m' and used for running DEFT-FUNNEL
% on a collection of test problems. It sets the objective function of the
% test problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( size( varargin, 2 ) > 0 )
   nprob = varargin{ 1 };
else
   nprob = 1;
end

if ( nprob == 1 )
    
   % Problem HS6
   f = problem_hs6_obj( x );
   
elseif ( nprob == 2 )
    
   % Problem HS7
   f = problem_hs7_obj( x );
   
elseif ( nprob == 3 )
    
   % Problem HS9
   f = problem_hs9_obj( x );
   
elseif ( nprob == 4 )
    
   % Problem HS10
   f = problem_hs10_obj( x );
   
elseif ( nprob == 5 )
    
   % Problem HS12
   f = problem_hs12_obj( x );
   
elseif ( nprob == 6 )
    
   % Problem HS13
   f = problem_hs13_obj( x );
   
elseif ( nprob == 7 ) 
   
   % Problem HS21
   f = problem_hs21_obj( x );

elseif ( nprob == 8 )
    
   % Problem HS23
   f = problem_hs23_obj( x );
   
elseif ( nprob == 9 )
    
   % Problem BT1
   f = problem_bt1_obj( x );
   
elseif ( nprob == 10 )
    
   % HANDBOOK - QUADRATICALLY CONSTRAINED TEST PROB 3
   f = problem_handbook_quadcons_pb3_obj( x );
   
end

