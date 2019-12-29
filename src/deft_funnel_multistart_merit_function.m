function [merit_fun_S, cons_viol_S] =                                       ...
    deft_funnel_multistart_merit_function(S, fS, cS, hS, nb_cons_c,         ...
    nb_cons_h, ls, us)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Computes the merit function value and the constraint vioation of every 
% point in S.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbsamples   = size(S, 2);
penalty     = 50;
merit_fun_S = zeros(nbsamples, 1);
cons_viol_S = zeros(nbsamples, 1);

for i = 1:nbsamples
    
    % Black-box constraints c
    for j=1:nb_cons_c
        if ~isinf(us(j))
            cons_viol_S(i) = cons_viol_S(i) + penalty * max(0, cS(j,i)-us(j));
        end
        if ~isinf(ls(j))
            cons_viol_S(i) = cons_viol_S(i) + penalty * max(0, ls(j)-cS(j,i));
        end
    end
    
    % White-box constraints h
    for j=1:nb_cons_h
        if ~isinf(us(nb_cons_c+j))
            cons_viol_S(i) = cons_viol_S(i) + penalty * max(0, hS(j,i)-us(nb_cons_c+j));
        end
        if ~isinf(ls(nb_cons_c+j))
            cons_viol_S(i) = cons_viol_S(i) + penalty * max(0, ls(nb_cons_c+j)-hS(j,i));
        end
    end
    merit_fun_S(i) = fS(i) + cons_viol_S(i);
    
end

end % end of deft_funnel_merit_function
