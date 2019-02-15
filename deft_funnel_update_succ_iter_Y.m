function [ sampleSet, setting, xstatus, pos ] =                             ...
             deft_funnel_update_succ_iter_Y( sampleSet, iterate_plus,       ...
             setting, modelSize, xstatus, const, succ )
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Desc: Tries to include the successful trial point into the interpolation set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Augment interpolation set if not fully quadratic yet
if ( setting.cur_degree < modelSize.pquad || ( setting.whichmodel == 3 &&   ...
                             setting.cur_degree < 2*modelSize.pquad ) )
    
    [ sampleSet, pY ] = deft_funnel_augment_Y( sampleSet, iterate_plus.x, setting );
    setting.cur_degree = pY;
    pos = pY;
    
else

    % Include xplus in the interpolation set, by replacing
    % another point if the model is already fully quadratic.
    [ sampleSet, pos ] = deft_funnel_include_in_Y( sampleSet,               ...
        iterate_plus.x, setting, [1:setting.cur_degree], setting.Lambda_XN, ...
        setting.criterion_S, succ );
    
    if ( pos > 0 )
        xstatus( sampleSet.ind_Y( pos ) )  = const.unused;
    end
    
end

if ( pos > 0 )
    xstatus( sampleSet.nbPoints ) = const.inY;
    sampleSet.ind_Y( pos )        = sampleSet.nbPoints;
    sampleSet.fY( pos )           = iterate_plus.feval;
    sampleSet.cY( :, pos )        = iterate_plus.ceval;
    if ( setting.verbose >= 2 )
        disp( [' replacement/inclusion of interpolation point at position', ...
            int2str( pos ), ' in Y successful'] )
        disp( ' ' )
    end
end
    
end % end of deft_funnel_update_succ_iter_Y