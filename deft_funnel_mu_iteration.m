function [ sampleSet, evaluations, Delta_f, Delta_c, xstatus, sstatus,      ...
           it_type, poised_model ] = deft_funnel_mu_iteration( f, c,        ...
             sampleSet, iterate, setting, evaluations, indicators, const,   ...
             xstatus, sstatus, exit_nstep, Delta_f, Delta_c, it_type, poised_model )
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Computations for the mu-iterations with possible change in the
% interpolation set in order to improve the surrogate models.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If the error in the models is not close to zero, improve the
% models by reducing the radius of the interpolation set.
if ( sampleSet.errg > setting.epsilon*1.0e-3 )

    oldErrg = sampleSet.errg;
    oldRadius = sampleSet.Y_radius;
    
    radius = 0.3*sampleSet.Y_radius;
    
    [ sampleSet, replaced ] = deft_funnel_repair_Y( sampleSet,         ...
         iterate, setting, radius);
     
    poised_model = 1;
    
    sample = deft_funnel_create_sample( [], [], iterate );
  
    % Compute the corresponding function values.
    for i = 1:length( replaced )
    
        j = replaced( i );

        % Set index of new point and update status of the old point
        xstatus( sampleSet.ind_Y( j ) ) = const.unused;
        sampleSet.ind_Y( j ) = sampleSet.nbPoints+1;
        sample.x = sampleSet.Y( :, j );

        % Update X and evaluate function at ynew
        [ sampleSet, sample, evaluations, xstatus, sstatus ] =          ...
        deft_funnel_augmX_evalfc( f, c, sample, sampleSet, setting,     ...
        evaluations, xstatus, const.inY, sstatus );

        % Update Y function values
        sampleSet.fY( j )    = sampleSet.fX( sampleSet.nbPoints );
        sampleSet.cY( :, j ) = sampleSet.cX( :, sampleSet.nbPoints );

    end

    if ( setting.verbose >= 2 )
        if ( oldRadius > sampleSet.Y_radius )
            disp( ' Y radius reduced to improve poisedness. ' )
            disp( ' ' )
        end
    end
    
    if ( sampleSet.errg < oldErrg )
        
        % If poisedness was improved, allow the trust region to increase in the
        % same order of the improvement
        it_type = strcat( it_type, ': improving poisedness' );
        Delta_f = Delta_f * (1 + ( oldErrg - sampleSet.errg));
        Delta_c = Delta_c * (1 + ( oldErrg - sampleSet.errg));
        
    elseif ( sampleSet.errg > setting.epsilon_i )
        Delta_f = Delta_f * 3*setting.gamma3;
        Delta_c = Delta_c * 3*setting.gamma3;
    end

elseif ( indicators.norm_c_s > setting.epsilon && exit_nstep == 1 )
    
    % Pathological case:
    % The current point is infeasible, but BLLs was not able to
    % compute a good normal step although the error in the interpolation
    % model is small. Try reducing the trust region radius.
    Delta_c = setting.gamma2 * Delta_c;
    
end
    
end % end of deft_funnel_mu_iteration