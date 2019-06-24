function [ sampleSet, iterate, setting, evaluations, xstatus, sstatus,      ...
           Delta_f, Delta_c, Delta_f_counter, vmax, rho, pos, it_type,      ...
           poised_model ] =                                                 ...
             deft_funnel_sec_order_correction( f, c, sampleSet, iterate,    ...
             iterate_plus, setting, models, indicators, evaluations,        ...
             d, delta_f, vmax, pos, xstatus, sstatus, const, modelSize,     ...
             Delta, Delta_f, Delta_c, Deltamax, Delta_f_counter, rho,       ...
             succ, it_type, poised_model )
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Desc: Calculates a second-order correction for the normal step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = iterate.xdim;
m = iterate.sdim;

nmax = min( Delta_c, setting.kappa_n * norm(iterate_plus.c_s) );
lb( 1:n+m )   = -nmax;
ub( 1:n+m )   = nmax;
dlx           = setting.lx( iterate.indfree ) - iterate_plus.x;
dlx           = dlx';
dux           = setting.ux( iterate.indfree ) - iterate_plus.x;
dux           = dux';
lb( 1:n )     = max( lb( 1:n ), dlx );
ub( 1:n )     = min( ub( 1:n ), dux );

dls           = setting.ls - iterate_plus.s;
dls           = dls';
dus           = setting.us - iterate_plus.s;
dus           = dus';
lb( n+1:n+m ) = max( lb( n+1:n+m ), dls );
ub( n+1:n+m ) = min( ub( n+1:n+m ), dus );

full_n  = length( iterate.xfix );
I = eye( full_n );

% Compute the second-order correction
sol      = deft_funnel_blls_exp( models.derivatives.J_s, -iterate_plus.c_s, lb', ub' );
n_xsoc   = sol( 1 : n );
n_ssoc   = sol( n+1:n+m );
norm_soc = norm( d + sol );

% Initialize the 'second-order-correction' point
iterate_soc = iterate_plus;

if ( norm_soc <= Delta )

    iterate_soc.x = iterate_plus.x + n_xsoc ;
    iterate_soc.s = iterate_plus.s + n_ssoc ;
    
    % Error message in case of evaluation failure
    error_msg = [ ' Possible causes:\n', ...
              ' (1) Hidden constraint: your black-box function might not',  ...
              ' be able to be evaluated at that point.\n',                  ...
              ' (2) Function handle: check if the function handle',         ...
              ' that you have passed as input is correct.\n',               ...
              ' (3) Black-box function: check if your black-box',           ...
              ' function does not contain any internal errors.\n',          ...
              ' Setting output(s) to 1.0e+10. This might affect the',       ...
              ' final solution.\n\n' ];

    if ( iterate.nfix > 0 )
        xsocfull = I( :, iterate_soc.indfix )*iterate_soc.xfix( iterate_soc.indfix ) + ...
            I( :, iterate_soc.indfree )*iterate_soc.x;

        if ( setting.scaleX )
            xsocfull = xsocfull ./ setting.scalefacX;
        end
        
        if ( strcmp( c, 'combined' ) )

            try
                output = f(xsocfull);
                iterate_soc.feval = output(1);
                iterate_soc.ceval = output(2:iterate.sdim+1);
            catch
                disp(' ')
                disp( ' Error: evaluation of the black box FAILED at the point');
                xsocfull
                fprintf(error_msg);
                iterate_soc.feval = 1.0e+10;
                iterate_soc.ceval = 1.0e+10 * ones( iterate.sdim, 1 );
            end

        else
            
            try
                iterate_soc.feval = f(xsocfull);
            catch
                disp(' ')
                disp(' Error: evaluation of objective function FAILED at the point' );
                xsocfull
                fprintf(error_msg);
                iterate_soc.feval = 1.0e+10;
            end

            try
                iterate_soc.ceval = c(xsocfull);
            catch
                disp(' ')
                disp(' Error: evaluation of constraint function FAILED at the point' );
                xsocfull
                fprintf(error_msg);
                iterate_soc.ceval = 1.0e+10 * ones( iterate.sdim, 1 );
            end

        end

    else

        if ( setting.scaleX )
            iterate_soc.x = iterate_soc.x ./ setting.scalefacX;
        end
        
        if ( strcmp( c, 'combined' ) )

            try
                output = f(iterate_soc.x);
                iterate_soc.feval = output(1);
                iterate_soc.ceval = output(2:iterate.sdim+1);
            catch
                disp(' ')
                disp( ' Error: evaluation of the black box FAILED at the point');
                iterate_soc.x
                fprintf(error_msg);
                iterate_soc.feval = 1.0e+10;
                iterate_soc.ceval = 1.0e+10 * ones( iterate.sdim, 1 );
            end

        else
            
            try
                iterate_soc.feval = f(iterate_soc.x);
            catch
                disp(' ')
                disp(' Error: evaluation of objective function FAILED at the point' );
                iterate_soc.x
                fprintf(error_msg);
                iterate_soc.feval = 1.0e+10;
            end

            try
                iterate_soc.ceval = c(iterate_soc.x);
            catch
                disp(' ')
                disp(' Error: evaluation of constraint function FAILED at the point' );
                iterate_soc.x
                fprintf(error_msg);
                iterate_soc.ceval = 1.0e+10 * ones( iterate.sdim, 1 );
            end

        end
        
    end

    %  Update evaluation counter
    evaluations.nfeval = evaluations.nfeval + 1;
    evaluations.nceval = evaluations.nceval + 1;
    if ( size( iterate_soc.ceval, 2 ) > 1 )
        iterate_soc.ceval = iterate_soc.ceval';
    end
    iterate_soc.c_s = iterate_soc.ceval - iterate_soc.s;
    vsoc = 0.5 * ( iterate_soc.c_s.' * iterate_soc.c_s );

    rho = ( iterate.feval - iterate_soc.feval + setting.rho_eps ) /         ...
          ( delta_f + setting.rho_eps );

    if ( rho >= setting.eta1 && vsoc <= vmax )

        Delta_f_counter = 0;

        % Set index of new point
        sampleSet.nbPoints = sampleSet.nbPoints + 1;

        %  Set initial status of the new point
        xstatus( sampleSet.nbPoints ) = 0;
        %  The point is contained in current subspace
        sstatus( sampleSet.nbPoints ) = 1;

        % Augment X with the new point

        if ( iterate_soc.nfix > 0 )
            sampleSet.X( :, sampleSet.nbPoints ) = xsocfull;
        else
            sampleSet.X( :, sampleSet.nbPoints ) = iterate_soc.x;
        end

        % Augment fX with new function value
        if ( max( size( iterate_soc.feval ) ) > 1 )
            sampleSet.fX( sampleSet.nbPoints ) = min( setting.fxmax,        ...
                norm( iterate_soc.feval )^2);
        else
            sampleSet.fX( sampleSet.nbPoints ) = min( setting.fxmax,        ...
                real( iterate_soc.feval ) );
        end

        % Augment cX with new function value
        length_c = max( size( iterate_soc.ceval ) );
        if ( length_c > 1 )
            for i=1:length_c
                sampleSet.cX( i, sampleSet.nbPoints ) = min( setting.fxmax, ...
                    real( iterate_soc.ceval( i ) ) );
            end
        else
            sampleSet.cX( sampleSet.nbPoints ) = min( setting.fxmax,        ...
                real( iterate_soc.ceval ) );
        end

        if ( setting.verbose >=2 )
            disp( ' Successful f-iteration with 2nd-order correction :' )
            disp( ' ' )
        end

        it_type = 'f-succ(2nd order correction used)';
        
        [ sampleSet, iterate, setting, pos, xstatus, Delta_f,               ...
          Delta_c, Delta_f_counter, vmax, poised_model ] =                  ...
            deft_funnel_succ_iteration( 1, sampleSet, iterate,              ...
            iterate_soc, setting, indicators, modelSize, Delta_f,           ...
            Delta_c, Delta_f_counter, Deltamax, vsoc, vmax,                 ...
            xstatus, const, succ, norm_soc, norm_soc, rho, poised_model);
    end
end
    
end % end of deft_funnel_sec_order_correction