function [ best_sol, best_fval, best_indicators, total_eval,                ...
           nb_local_searches, fX ] = deft_funnel_multistart( objfun, cons,  ...
           n, nbcons, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative-Free Trust FUNNEL (DEFT-FUNNEL) for Global Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file implements the Multi-level Single Linkage (MLSL) method, which is
% a clustering-based multistart technique for global optimization. The
% MLSL is used to select the starting points of the local searches done by
% the SQP trust-region-based algorithm called in 'deft_funnel.m'.
%
% Output: this code returns the best feasible local minimum found under a
% predefined total number of simulations.
%
% Please read the README file before embarking on this journey.
% Last update: 11/02/2019.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
X = [];                     % set of stationary points
fX = [];                    % set of optimal values found
S = [];                     % set of sample points
started_points = [];        % indicates which points in S have been used
fS = [];                    % objective function values of points in S
cS = [];                    % constraint values of points in S
mS = [];                    % merit function values of points in S
best_sol = NaN*ones(n,1);   % best solution from all local searches
best_fval = Inf;            % best feasible objective function value
best_indicators = [];       % indicators of the best solution found
total_eval = 0;             % total number of evaluations
f_global_optimum = NaN;     % known objective function value of the global optimum
nb_global_found = 0;        % number of times the global optimum was found
nb_local_searches = 0;      % number of local searches run
maxeval = 5000*n;           % maximum number of simulations (budget)
maxeval_ls = maxeval*0.7;   % maximum number of simulations per local search
tol_feas = 1.0e-4;          % tolerance for feasibility
tol_conv = 1.0e-4;          % tolerance for declaring convergence based on the 
                            % distance to the known global optimum (not
                            % used if the global optimum is not known a
                            % priori)

% Set the MLSL constants
nb_of_new_samples = min(20*n, floor(0.1*maxeval));
sigma = 4;
kappa = 1;                  % percentage of the best samples to consider
maxit = 20;                 % number of MLSL iterations
alpha = pi^(-0.5) * ( gamma(1+n/2) * sigma )^(1/n); % MLSL formula

% Setting default bounds
lx(1:n) = -Inf;             % Default 'x' lower bounds
lx=lx';      
ux(1:n) = Inf;              % Default 'x' upper bounds
ux=ux';      
ls(1:nbcons) = 0.0;         % Default 's' lower bounds
ls=ls';      
us(1:nbcons) = 0.0;         % Default 's' upper bounds
us=us';  

% Check the parity of the variable argument list.
noptargs = size( varargin, 2 );    % The number of optional arguments
if mod( noptargs, 2 ) > 0
   msg   = [' DEFT-FUNNEL error: the number of variable arguments',         ...
            ' must be even! Default parameters used.' ];
   disp( msg );
   return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%      PROCESS ARGUMENT LIST      %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:2:noptargs

    % Lower bounds on the x variables
    if ( strcmp( varargin{ i }, 'lxbounds') )
        if ( isnumeric( varargin{ i + 1 } ) )
            if ( isempty( varargin{ i + 1 } ) == 0 && ...
                 length( varargin{ i + 1 } ) == n )
                lx = varargin{ i + 1 };
                if ( size( lx, 1 ) == 1 )
                    lx=lx';
                end
            else
                msg = [ ' DEFT-FUNNEL multistart error: lower x bounds',   	...
                        ' empty or not of length n! No bounds used.' ];
                disp( msg )
            end
        else
            msg = [ ' DEFT-FUNNEL multistart error: wrong type of input',  	... 
                    ' for lower bounds! No bounds used.'];
            disp( msg )
        end
        
    % Upper bounds on the x variables
    elseif ( strcmp( varargin{ i }, 'uxbounds' ) )
        if ( isnumeric( varargin{ i + 1 } ) )
            if ( isempty( varargin{ i + 1 } ) == 0 && ...
                 length( varargin{ i + 1 } ) == n )
                ux = varargin{ i + 1 };
                if ( size( ux, 1 ) == 1 )
                    ux = ux';
                end
            else
                msg = [ ' DEFT-FUNNEL multistart error: upper x bounds', 	...
                        ' empty or not of length n! No bounds used.' ];
                disp( msg )
            end
        else
            msg = [ ' DEFT-FUNNEL multistart error: wrong type of input', 	... 
                    ' for upper bounds! No bounds used.' ];
            disp( msg )
        end
        
    % Lower bounds on the s variables
    elseif ( strcmp( varargin{ i }, 'lsbounds' ) )
        if ( isnumeric( varargin{ i + 1 } ) )
            if ( isempty( varargin{ i + 1 } ) == 0 && ...
                 length( varargin{ i + 1 } ) == nbcons )
                ls = varargin{ i + 1 };
                if ( size( ls, 1 ) == 1 )
                    ls = ls';
                end
            else
                msg = [ ' DEFT-FUNNEL multistart error: lower s bounds',    ...
                        '  empty or not of length n! No bounds used.' ];
                disp( msg )
            end
        else
            msg = [ ' DEFT-FUNNEL multistart error: wrong type of input',   ... 
                    ' for lower bounds! No bounds used.' ];
            disp( msg )
        end
        
    % Upper bounds on the s variables
    elseif ( strcmp( varargin{ i }, 'usbounds' ) )
        if ( isnumeric( varargin{ i + 1 } ) )
            if ( isempty( varargin{ i + 1 } ) == 0 && ...
                 length( varargin{ i + 1 } ) == nbcons )
                us = varargin{ i + 1 };
                if ( size( us, 1 ) == 1 )
                    us = us';
                end
            else
                msg = [ ' DEFT-FUNNEL multistart error: upper s bounds',    ...
                        ' empty or not of length n! No bounds used.' ];
                disp( msg )
            end
        else
            msg = [ ' DEFT-FUNNEL multistart error: wrong type of input',   ... 
                    ' for upper bounds! No bounds used.' ];
            disp( msg )
        end
        
    % Maximum number of simulations
    elseif ( strcmp( varargin{ i }, 'maxeval' ) )
        if ( isnumeric( varargin{ i + 1 } ) )
            maxeval = varargin{ i + 1 };
        else
            msg = [ ' DEFT-FUNNEL error: wrong type of input for',          ... 
                    ' multistart entry!' ];
            disp( msg )
        end
        
    % Maximum number of simulations per local search
    elseif ( strcmp( varargin{ i }, 'maxeval_ls' ) )
        if ( isnumeric( varargin{ i + 1 } ) )
            maxeval_ls = varargin{ i + 1 };
        else
            msg = [ ' DEFT-FUNNEL error: wrong type of input for',          ... 
                    ' multistart entry!' ];
            disp( msg )
        end
        
    % Ojective function value of the global optimum known by the user
    elseif ( strcmp( varargin{ i }, 'f_global_optimum' ) )
        if ( isnumeric( varargin{ i + 1 } ) )
            maxeval_ls = varargin{ i + 1 };
        else
            msg = [ ' DEFT-FUNNEL error: wrong type of input for',          ... 
                    ' multistart entry!' ];
            disp( msg )
        end
        
    else
        msg = [ ' DEFT-FUNNEL multistart warning: undefined keyword',       ...
                varargin{ i }, '! Ignoring.' ];
        disp( msg )
    end
end

% Calculate a rough approximation of the volume of the feasible region.
% If the feasible region is not bounded by a box, set the 'x' lower/upper
% bounds to -10/10 when necessary
volume_factor = 1;
for i = 1:n
    a = lx(i);
    b = ux(i);
    if isinf(a)
        a = -10;
    end
    if isinf(b)
        b = max(0,a) + 10;
    end
    volume_factor = volume_factor*(b - a);
end
alpha = alpha * volume_factor^(1.0/n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%      MULTISTART SEARCH      %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:maxit
    
    if ( total_eval >= maxeval )
        disp(' TERMINATED: maximum budget reached')
        disp(' ')
        break
    end
    
    if ( total_eval + nb_of_new_samples > maxeval )
        total_eval
        maxeval
        disp(' TERMINATED: the remaining budget is not enough to continue')
        disp(' ')
        break
    end
    
    % Sample new points and add them to S
    Snew = deft_funnel_multistart_sampling('uniform', n, nb_of_new_samples, lx, ux);
    nbsamples_old = size( S, 2);
    S = [ S Snew ];
    started_points = [ started_points zeros(1, nb_of_new_samples) ];
    nbsamples = size( S, 2);
    
    % Calculate their objective, constraint and merit function values
    for i = 1:size( Snew, 2 )
       fS( nbsamples_old + i ) = objfun( Snew(:, i) );
       cS( :, nbsamples_old + i ) = cons( Snew(:, i) );
    end
    total_eval = total_eval + nb_of_new_samples;
    
    [ mSnew, cons_viol_Snew ] = deft_funnel_multistart_merit_function(Snew, ...
        fS(nbsamples_old + 1: end), cS(:,nbsamples_old + 1: end), ls, us);
    mS = [ mS mSnew ];
    
    [ Sreduced, mSreduced, asc_order ] = deft_funnel_multistart_best_samples( S, mS, kappa );
    
    % Update cluster radius
    rk = alpha * (log(k*nb_of_new_samples)/(k*nb_of_new_samples) )^(1/n);
        
    for j = 1:size( Sreduced, 2 )
        
        belong_cluster = deft_funnel_multistart_cluster( Sreduced, mSreduced, j, rk );
        
        % Check if there is a better point in a nearby cluster or if it
        % has already been used as a starting point
        if ( ~belong_cluster && ~started_points(asc_order(j)))
            
            started_points(asc_order(j)) = 1;
            x0 = Sreduced(:,j);
            total_eval = total_eval - 1;
            nb_local_searches = nb_local_searches + 1;
            
            % Start local search
            [x, fx, mu, indicators, evaluations, iterate, exit_algo ] =     ...
              deft_funnel( objfun, cons, x0, nbcons,                        ...
              'lxbounds', lx, 'uxbounds' , ux,                              ...
              'lsbounds', ls, 'usbounds' , us,                              ...
              'multistart_call', nb_local_searches,                         ...
              'maxeval' , min(maxeval_ls, maxeval - total_eval));
              
            if ( exit_algo == 0 ) % no errors found
                X = [X x];
                fX = [fX fx];
                if ( fx < best_fval && indicators.norm_c_s <= tol_feas )
                    best_fval = fx;
                    best_sol = x;
                    best_indicators = indicators;
                end
                
                if ( ~isnan(f_global_optimum) )
                    if ( abs(fx - f_global_optimum) < tol_conv )
                        nb_global_found = nb_global_found + 1;
                    end
                end
            end
           
            total_eval = total_eval + evaluations.nfeval;
            
        end
    end
end

end % end of deft_funnel_multistart

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sample points from the uniform distribution
function samples = deft_funnel_multistart_sampling(distribution, dimension, nbsamples, lx, ux)

if ( strcmp( distribution, 'uniform' ) )
    
    % Generate 'nbsamples' points between bounds 'lx' and 'ux'.
    % If the lower bound or the upper bound is not a real number, consider -10
    % and 10, respectively.
    rng('shuffle', 'twister');
    for i=1:dimension
        if ~isinf(lx(i))
            a = lx(i);
        else
            a = -10;
        end
        if ~isinf(ux(i))
            b = ux(i);
        else
            b = max(0,a) + 10;
        end
        samples(i,:) = (b-a).*rand(1,nbsamples) + a;
    end
end

end % end of deft_funnel_multistart_sampling