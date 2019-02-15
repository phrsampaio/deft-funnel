function A = deft_funnel_multistart_run_single_prob

% % Problem G6
% 
% A = [];
% X = [];
% 
% for i=1:30
%     
%     [ best_sol, best_fval, best_indicators, total_eval, nb_local_searches ] = deft_funnel_multistart( @problem_G6_obj, @problem_G6_cons, 2, 2, 'lsbounds', [ 0 0 ], 'usbounds', [ Inf Inf ], 'lxbounds', [ 13 0 ], 'uxbounds', [ 100 100 ]);
%     
%     if isempty(best_indicators)
%         norm_c_s = Inf;
%     else
%         norm_c_s = best_indicators.norm_c_s;
%     end
%     A = [ A; best_fval norm_c_s total_eval nb_local_searches];
%     X = [ X; best_sol' ];
% end
% 
% A = A(all(~isinf(A),2),:);
% 
% fileID = fopen('200_sim_test_prob_G6.txt','w');
% for i = 1:size(A,1)
%     fprintf(fileID,'%12.10f\t',A(i,:));
%     fprintf(fileID,'\n');
% end
% fclose(fileID);
% 
% fileID = fopen('100_sim_test_prob_G6_X.txt','w');
% for i = 1:size(X,1)
%     fprintf(fileID,'%4.8f\t',X(i,:));
%     fprintf(fileID,'\n');
% end
% fclose(fileID);

% Problem G8

% A = [];
% X = [];
% 
% for i=1:30
%     
%     [ best_sol, best_fval, best_indicators, total_eval, nb_local_searches ] = deft_funnel_multistart( @problem_G8_obj, @problem_G8_cons, 2, 2, 'lsbounds', [ -Inf -Inf ], 'usbounds', [ 0 0 ], 'lxbounds', [ 0 0 ], 'uxbounds', [ 10 10 ]);
%     
%     if isempty(best_indicators)
%         norm_c_s = Inf;
%     else
%         norm_c_s = best_indicators.norm_c_s;
%     end
%     A = [ A; best_fval norm_c_s total_eval nb_local_searches];
%     X = [ X; best_sol' ];
% end
% 
% A = A(all(~isinf(A),2),:);
% 
% fileID = fopen('50_sim_test_prob_G8.txt','w');
% for i = 1:size(A,1)
%     fprintf(fileID,'%12.10f\t',A(i,:));
%     fprintf(fileID,'\n');
% end
% fclose(fileID);
% 
% fileID = fopen('50_sim_test_prob_G8_X.txt','w');
% for i = 1:size(X,1)
%     fprintf(fileID,'%4.8f\t',X(i,:));
%     fprintf(fileID,'\n');
% end
% fclose(fileID);

% % Problem G9
% 
% A = [];
% X = [];
% 
% for i=1:1
%     
%     [ best_sol, best_fval, best_indicators, total_eval, nb_local_searches ] = deft_funnel_multistart( @problem_G9_obj, @problem_G9_cons, 7, 4, 'lsbounds', [ -Inf -Inf -Inf -Inf ], 'usbounds', [ 0 0 0 0 ], 'lxbounds', -10*ones(1,7), 'uxbounds', 10*ones(1,7));
%     
%     A = [ A; best_fval best_indicators.norm_c_s total_eval nb_local_searches];
%     X = [ X; best_sol' ];
% end
% 
% fileID = fopen('test_prob_G9.txt','w');
% for i = 1:size(A,1)
%     fprintf(fileID,'%12.10f\t',A(i,:));
%     fprintf(fileID,'\n');
% end
% fclose(fileID);
% 
% fileID = fopen('test_prob_G9_X.txt','w');
% for i = 1:size(X,1)
%     fprintf(fileID,'%4.8f\t',X(i,:));
%     fprintf(fileID,'\n');
% end
% fclose(fileID);

% Problem G11

% A = [];
% X = [];
% 
% for i=1:15
%     
%     [ best_sol, best_fval, best_indicators, total_eval, nb_local_searches ] = deft_funnel_multistart( @problem_G11_obj, @problem_G11_cons, 2, 1, 'lsbounds', 0, 'usbounds', 0, 'lxbounds', -1*ones(1,2), 'uxbounds', ones(1,2));
%     
%     A = [ A; best_fval best_indicators.norm_c_s total_eval nb_local_searches];
%     X = [ X; best_sol' ];
% end
% 
% fileID = fopen('test_prob_G11_2.txt','w');
% for i = 1:size(A,1)
%     fprintf(fileID,'%12.10f\t',A(i,:));
%     fprintf(fileID,'\n');
% end
% fclose(fileID);
% 
% fileID = fopen('test_prob_G11_X_2.txt','w');
% for i = 1:size(X,1)
%     fprintf(fileID,'%4.8f\t',X(i,:));
%     fprintf(fileID,'\n');
% end
% fclose(fileID);

% Problem Test problem 3

A = [];
X = [];

for i=1:30
    
    [ best_sol, best_fval, best_indicators, total_eval, nb_local_searches ] = deft_funnel_multistart( @handbook_quadcons_pb3_obj, @handbook_quadcons_pb3_cons, 6, 5, 'lsbounds', [ 4 4 -Inf -Inf 2], 'usbounds', [ Inf Inf 2 2 6 ], 'lxbounds', [ 0 0 1 0 1 0 ], 'uxbounds', [ Inf Inf 5 6 5 10 ]);
    
    if isempty(best_indicators)
        norm_c_s = Inf;
    else
        norm_c_s = best_indicators.norm_c_s;
    end
    A = [ A; best_fval norm_c_s total_eval nb_local_searches];
    X = [ X; best_sol' ];
end
A = A(all(~isinf(A),2),:);

fileID = fopen('25_sim_test_prob_handbook_quadcons_pb3.txt','w');
for i = 1:size(A,1)
    fprintf(fileID,'%12.10f\t',A(i,:));
    fprintf(fileID,'\n');
end
fclose(fileID);

fileID = fopen('25_sim_test_prob_handbook_quadcons_pb3_X.txt','w');
for i = 1:size(X,1)
    fprintf(fileID,'%4.8f\t',X(i,:));
    fprintf(fileID,'\n');
end
fclose(fileID);

% % Problem PrP
% 
% A = [];
% X = [];
% 
% for i=1:15
%     
%     [ best_sol, best_fval, best_indicators, total_eval, nb_local_searches ] = deft_funnel_multistart( @problem_PrP_obj, @problem_PrP_cons, 4, 3, 'lsbounds', [ -Inf -Inf -Inf ], 'usbounds', [ 0 0 0 ], 'lxbounds', [ 0 0 0 0 ], 'uxbounds', [ 1 1 50 240 ] );
%     
%     if isempty(best_indicators)
%         norm_c_s = Inf;
%     else
%         norm_c_s = best_indicators.norm_c_s;
%     end
%     A = [ A; best_fval norm_c_s total_eval nb_local_searches];
%     X = [ X; best_sol' ];
% end
% 
% A = A(all(~isinf(A),2),:);
% 
% fileID = fopen('test_prob_PrP_2.txt','w');
% for i = 1:size(A,1)
%     fprintf(fileID,'%12.10f\t',A(i,:));
%     fprintf(fileID,'\n');
% end
% fclose(fileID);
% 
% fileID = fopen('test_prob_PrP_X_2.txt','w');
% for i = 1:size(X,1)
%     fprintf(fileID,'%4.8f\t',X(i,:));
%     fprintf(fileID,'\n');
% end
% fclose(fileID);

end % end of deft_funnel_multistart_run_single_prob