function [ x, fx, mu, indicators, evaluations, iterate, exit_algo ] =       ...
    deft_funnel( f, c, x0, nbcons, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative-Free Trust FUNNEL (DEFT-FUNNEL) for Local Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Global Optimization, check the file 'deft_funnel_multistart.m'.
% Please read the README file before embarking on this journey.
% Last update: 15/05/2019.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check the type of mandatory inputs.
if ( ( ~isa( c, 'function_handle' ) && ~strcmp( c, 'combined' ) ) ||        ...
     ~isa( f, 'function_handle' ) || ~isnumeric( x0 ) || ~isnumeric( nbcons ) )
    
    if ( ~isa( c, 'function_handle' ) && ~strcmp( c, 'combined' ) )
        msg = [ ' DEFT-FUNNEL LOCALOPT error: the second argument is',      ...
                ' neither a function handle nor the string "combined"!',    ...
                ' Terminating.' ];
    end
    if ( ~isa( f, 'function_handle' ) )
        msg = [ ' DEFT-FUNNEL LOCALOPT error: the first argument is not a', ...
                ' function handle! Terminating.' ];
    end
    if ( ~isnumeric( x0 ) )
        msg = [ ' DEFT-FUNNEL LOCALOPT error: the starting point is not a', ...
                ' numeric vector! Terminating.' ];
    end
    if ( ~isnumeric( nbcons ) )
        msg = [ ' DEFT-FUNNEL LOCALOPT error: the number of constraints',   ...
                ' is not numeric! Terminating.' ];
    end
    disp( msg )
    return
end

% Reset the random generator 
rng('shuffle', 'twister');

% Make sure the starting point is a column vector
if ( size( x0, 1 ) == 1 && size( x0, 2 ) > 1 )
   x0 = x0';
end

% Set the dimension of the problem
n = length( x0 );

% Check for faulty input dimensions 
if( nbcons < 0 || n < 0 )
  msg = ' DEFT-FUNNEL LOCALOPT error: The problem dimensions are faulty.';
  disp( msg )
  return
end

% Check the parity of the variable argument list.
noptargs = size( varargin, 2 ); % The number of optional arguments
if mod( noptargs, 2 ) > 0
   msg   = [' DEFT-FUNNEL LOCALOPT error: the number of variable',          ...
            ' arguments must be even! Default parameters used.' ];
   disp( msg );
   return
end

% Set some constants
const.free         = 0;              % Variable free
const.fixed        = 1;              % Variable fixed
const.alwaysfixed  = 2;              % Variable always fixed
const.in           = 1;              % Point is contained in the current subspace
const.out          = 0;              % Point is not contained in the current subspace
const.unused       = 0;              % Point is not in the current interpolation set Y
const.inY          = 1;              % Point is in the current interpolation set Y

% Define initial values
nit                = 0;
nitold             = 0;
sampleSet.X        = [];             % Set of all points already evaluated
sampleSet.fX       = [];             % Objective function values of points in X
sampleSet.cX       = [];             % Constraint functions' values of points in X
sampleSet.Y        = [];             % Set of interpolation points currently used
sampleSet.fY       = [];             % Objective function values of points in Y
sampleSet.cY       = [];             % Constraint functions' values of points in Y
sampleSet.nbPoints = 0;              % Total number of evaluated points, i.e. |X|.
                                     % Updated in 'deft_funnel_augmX_evalfc' only
sampleSet.ind_Y    = [];             % Indices of points of the interpolation set Y in X
                                     % Updated after calling 'deft_funnel_augmX_evalfc'
sampleSet.i_xbest  = NaN;            % Indice of the current iterate (i.e., Y(1)) in X
xstatus            = [];             % Indicates whether the points are in Y or not (1 or 0)
sstatus            = [];             % Indicates whether the points are in the subspace or not (1 or 0)
vstatus            = zeros( n, 1 );  % Variable status (const.free or const.fixed)
sspace_save        = [];             % All subspaces already explored
xspace_save        = [];             % Indices of entering and exiting points corresponding  
                                     % to each subspace already explored

modelSize.pquad    = ( ( n + 1 )*( n + 2 ) ) / 2;  % Size of a fully quadratic model
modelSize.pdiag    = 2 * n + 1;                    % Size of a diagonal model
modelSize.plin     = n + 1;                        % Size of a linear model
evaluations.nfeval = 0;
evaluations.nceval = 0;

% Set initial default trust regions radii
Delta_f0           = 1;              % The initial TR radius for the model of f
Delta_c0           = 1;              % The initial TR radius for the model of c

% Set fixed parameters
setting            = deft_funnel_set_parameters( n, nbcons, ...
                     modelSize.plin, modelSize.plin, 'simplex'); 

% Set the 'iterate' structure
iterate = deft_funnel_create_sample( x0, nbcons, [] );
iterate.X = iterate.x'; % For tracking iterates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%      PROCESS ARGUMENT LIST      %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:2:noptargs

    % Lower bounds on the x variables
    if ( strcmp( varargin{ i }, 'lxbounds') )
        if ( isnumeric( varargin{ i + 1 } ) )
            if ( isempty( varargin{ i + 1 } ) == 0 && ...
                 length( varargin{ i + 1 } ) == n )
                setting.lx = varargin{ i + 1 };
                if ( size( setting.lx, 1 ) == 1 )
                    setting.lx=setting.lx';
                end
            else
                msg = [ ' DEFT-FUNNEL LOCALOPT warning: lower x bounds',    ...
                        ' empty or not of length n! No bounds used.' ];
                disp( msg )
            end
        else
            msg = [ ' DEFT-FUNNEL LOCALOPT warning: wrong type of input',   ... 
                    ' for lower x bounds! No bounds used.'];
            disp( msg )
        end
        
    % Upper bounds on the x variables
    elseif ( strcmp( varargin{ i }, 'uxbounds' ) )
        if ( isnumeric( varargin{ i + 1 } ) )
            if ( isempty( varargin{ i + 1 } ) == 0 && ...
                 length( varargin{ i + 1 } ) == n )
                setting.ux = varargin{ i + 1 };
                if ( size( setting.ux, 1 ) == 1 )
                    setting.ux = setting.ux';
                end
            else
                msg = [ ' DEFT-FUNNEL LOCALOPT warning: upper x bounds',    ...
                        ' empty or not of length n! No bounds used.' ];
                disp( msg )
            end
        else
            msg = [ ' DEFT-FUNNEL LOCALOPT warning: wrong type of input',   ... 
                    ' for upper x bounds! No bounds used.' ];
            disp( msg )
        end
        
    % Lower bounds on the s variables
    elseif ( strcmp( varargin{ i }, 'lsbounds' ) )
        if ( isnumeric( varargin{ i + 1 } ) )
            if ( isempty( varargin{ i + 1 } ) == 0 && ...
                 length( varargin{ i + 1 } ) == nbcons )
                setting.ls = varargin{ i + 1 };
                if ( size( setting.ls, 1 ) == 1 )
                    setting.ls = setting.ls';
                end
            else
                msg = [ ' DEFT-FUNNEL LOCALOPT warning: lower s bounds',    ...
                        ' empty or not of length "nbcons"! Setting to',     ...
                        ' zero by default.' ];
                disp( msg )
            end
        else
            msg = [ ' DEFT-FUNNEL LOCALOPT warning: wrong type of input',   ... 
                    ' for lower s bounds! Setting to zero by default.' ];
            disp( msg )
        end
        
    % Upper bounds on the s variables
    elseif ( strcmp( varargin{ i }, 'usbounds' ) )
        if ( isnumeric( varargin{ i + 1 } ) )
            if ( isempty( varargin{ i + 1 } ) == 0 && ...
                 length( varargin{ i + 1 } ) == nbcons )
                setting.us = varargin{ i + 1 };
                if ( size( setting.us, 1 ) == 1 )
                    setting.us = setting.us';
                end
            else
                msg = [ ' DEFT-FUNNEL LOCALOPT warning: upper s bounds',    ...
                        ' empty or not of length "nbcons"! Set to zero by', ...
                        ' default.' ];
                disp( msg )
            end
        else
            msg = [ ' DEFT-FUNNEL LOCALOPT warning: wrong type of input',   ... 
                    ' for upper s bounds! Set to zero by default.' ];
            disp( msg )
        end
        
    % Check if multistart method is being used
    elseif ( strcmp( varargin{ i }, 'multistart_call' ) )
        if ( isnumeric( varargin{ i + 1 } ) )
            setting.multistart_call = varargin{ i + 1 };
        else
            msg = [ ' DEFT-FUNNEL LOCALOPT warning: wrong type of input',   ... 
                    ' for "multistart_call" entry!' ];
            disp( msg )
        end
    
    % Maximum number of simulations
    elseif ( strcmp( varargin{ i }, 'maxeval' ) )
        if ( isnumeric( varargin{ i + 1 } ) )
            setting.maxeval = varargin{ i + 1 };
        else
            msg = [ ' DEFT-FUNNEL LOCALOPT warning: wrong type of input',   ... 
                    ' for maximum number of simulations. Default value',    ...
                    ' will be used' ];
            disp( msg )
        end
        
    else
        msg = [ ' DEFT-FUNNEL LOCALOPT warning: undefined keyword',         ...
                varargin{ i }, '! Ignoring input value.' ];
        disp( msg )
    end
end

% Compute the maximal TR radius
Delta_f = Delta_f0;
Delta_c = Delta_c0;
Delta = min( Delta_f, Delta_c );
Deltamax = setting.factor_Dmax * Delta;

if ( setting.multistart_call <= 1 )
    deft_funnel_print_banner( setting );
else
   disp(' ')
   disp([' Local search number: ', int2str(setting.multistart_call)])
end

% Check the bounds and correct intial Delta if there is no sufficient space 
% between the bounds. Shift x0 if it is out of the box.
disp(' ')
for j=1:n
    
    % Check lower and upper bounds
    if setting.lx( j ) > setting.ux( j )
        msg  = [ ' DEFT-FUNNEL LOCALOPT warning: Lower bound of component ', ...
                 int2str( j ), ' exceeds upper bound !!' ];
        disp( msg )
        exit_algo = 1;
        return
    end
    
    % Check difference between bounds
    temp( j ) = setting.ux( j ) - setting.lx( j );
    
    if ( temp( j ) < Delta + Delta )
        if ( temp( j ) == 0.0 )
            iterate.nfix      = iterate.nfix + 1;
            iterate.indfix    = [ iterate.indfix j ];
            iterate.xfix( j ) = setting.lx( j );
            vstatus( j )      = const.alwaysfixed;
        else
            Delta = 0.5 * temp( j );
            [ Delta_f, Delta_c ] = deal( Delta );
            
            disp( [ ' Diff. between lower and upper bound of component ',   ...
                    int2str( j ), ' is less than 2*Delta0! New Delta0',     ...
                    ' set to ', num2str( Delta ) ] );
        end
    end
    
    % Move the starting point inside the bounds if necessary
    if ( iterate.x( j ) < setting.lx( j ) )
        iterate.x( j ) = setting.lx( j );
    elseif ( iterate.x( j ) > setting.ux( j ) )
        iterate.x( j ) = setting.ux( j );
    end
    
end

% Scale initial point and bounds, if user-defined
if ( setting.scaleX )
    for i = 1:n
        if ( setting.scalefacX( i ) > 0 )
            iterate.x(i)  = iterate.x( i ) * setting.scalefacX( i );
            setting.lx(i) = setting.lx( i ) * setting.scalefacX( i );
            setting.ux(i) = setting.ux( i ) * setting.scalefacX( i );
        else
            setting.scalefacX( i ) = 1;
        end
    end
end

% Reset constants if some variables are fixed
if ( iterate.nfix > 0 )
    nfree = iterate.fulldim - iterate.nfix;
    if ( nfree <= 0 )
        msg = ' DEFT-FUNNEL LOCAL error: No free variables. Please, enlarge search space.';
        disp( msg );
        exit_algo = 1;
        return
    end
    iterate.indfree = setdiff( 1:iterate.fulldim, iterate.indfix );
    iterate.x       = iterate.x( iterate.indfree );
    iterate.xdim    = nfree;
    if ( setting.cur_degree == modelSize.plin )
        setting.cur_degree = nfree + 1;
    elseif ( setting.cur_degree == modelSize.pdiag )
        setting.cur_degree = 2*nfree+1;
    elseif ( setting.cur_degree == modelSize.pquad )
        setting.cur_degree = ( ( nfree + 1 ) * ( nfree + 2 ) ) / 2;
    end
    if ( setting.rep_degree == modelSize.plin )
        setting.rep_degree = nfree+1;
    elseif ( setting.rep_degree == modelSize.pdiag )
        setting.rep_degree = 2*nfree+1;
    elseif ( setting.rep_degree == modelSize.pquad )
        setting.rep_degree = ( ( nfree + 1 ) * ( nfree + 2 ) ) / 2;
    end
    modelSize.plin  = nfree + 1;
    modelSize.pdiag = 2 * nfree + 1;
    modelSize.pquad = ( ( nfree + 1 ) * ( nfree + 2 ) ) / 2;
else
    iterate.indfree = 1:iterate.fulldim;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%      BUILDING THE INTIAL POISED SAMPLE SET      %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ sampleSet, iterate, evaluations, xstatus, sstatus, poised_model,          ...
  msg ] = deft_funnel_build_initial_sample_set( f, c, sampleSet,            ...
    iterate, setting, evaluations, xstatus, sstatus, const, modelSize,      ...
    Delta );

if ( strcmp( msg(2:6), 'Error' ) )
    disp(' ')
    disp(' DEFT-FUNNEL error: Construction of initial sample set unsuccessful.' )
    x = iterate.x;
    fx = iterate.feval;
    mu = NaN;
    indicators = NaN;
    exit_algo = 1;
    return
end

% Compute the maximal objective function value
setting.fxmax = min( setting.fmax, setting.factor_fmax * max( abs( iterate.feval ), 1 ) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%      INITIALIZATION OF THE ALGORITHM      %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set initial values for the slacks
iterate.s = zeros( iterate.sdim , 1 );
minThreshold = -1.0e+10;
maxThreshold = 1.0e+10;
for i=1:nbcons
    if ( isfinite( setting.ls( i ) ) && isfinite( setting.us( i ) ) )
        if ( setting.ls( i ) > minThreshold && setting.us( i ) < maxThreshold )
            if ( iterate.ceval( i ) > setting.ls( i ) && iterate.ceval( i ) < setting.us( i ) )
                iterate.s( i ) = iterate.ceval( i );
            else
                iterate.s( i ) = (setting.ls( i ) + setting.us( i ))/2;
            end
        elseif ( setting.ls( i ) > minThreshold )
            if ( iterate.ceval( i ) > setting.ls( i ) && iterate.ceval( i ) < setting.us( i ) )
                iterate.s( i ) = iterate.ceval( i );
            else
                iterate.s( i ) = setting.ls( i );
            end
        else
            if ( iterate.ceval( i ) > setting.ls( i ) && iterate.ceval( i ) < setting.us( i ) )
                iterate.s( i ) = iterate.ceval( i );
            else
                iterate.s( i ) = setting.us( i );
            end
        end
    elseif ( isfinite( setting.ls( i ) ) )
        if ( setting.ls( i ) > minThreshold )
            if ( iterate.ceval( i ) > setting.ls( i ) && iterate.ceval( i ) < setting.us( i ) )
                iterate.s( i ) = iterate.ceval( i );
            else
                iterate.s( i ) = setting.ls( i );
            end
        else
            if ( iterate.ceval( i ) > setting.ls( i ) && iterate.ceval( i ) < setting.us( i ) )
                iterate.s( i ) = iterate.ceval( i );
            else
                iterate.s( i ) = setting.us( i );
            end
        end
    else
        if ( iterate.ceval( i ) > setting.ls( i ) && iterate.ceval( i ) < setting.us( i ) )
            iterate.s( i ) = iterate.ceval( i );
        else
            iterate.s( i ) = setting.us( i );
        end
    end
end

if ( norm( iterate.ceval - iterate.s ) > 1.0e+3 )
    iterate.s = iterate.ceval;
end

infeasls = find( iterate.ceval < setting.ls );
if ( infeasls > 0 )
    iterate.s(infeasls) = setting.ls(infeasls);
end

infeasus = find( iterate.ceval > setting.us );
if ( infeasus > 0 )
    iterate.s(infeasus) = setting.us(infeasus);
end

% Compute Lagrange multipliers estimates and the gradient of the 
% Lagrangian function at the iterate. Build the models and get
% optimality and feasibility indicators.
[ models, iterate, indicators ] = ...
   deft_funnel_update_iteration( sampleSet, iterate, setting );

% Define the initial funnel radius
vmax = max( setting.kappa_ca, setting.kappa_cr * indicators.v );

% Define the initial convergence threshold
setting.epsilon_i = setting.epsilon0;

% Print the information obtained after computing the first model
deft_funnel_print_head_info( setting );
deft_funnel_printout( nit, evaluations, iterate, setting, ...
    indicators, vmax, 0.0, 0.0, Delta, 0.0, '', '', sampleSet.errg );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%      CALL MAIN ROUTINE      %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[  exit_algo, nit, sampleSet, iterate, indicators, evaluations ] =          ...
    deft_funnel_main( f, c, nit, nitold, sampleSet, iterate, setting,       ...
    indicators, evaluations, models, xstatus, sstatus, vstatus,             ...
    sspace_save, xspace_save, const, modelSize, Delta, Delta_f, Delta_c,    ...
    Deltamax, vmax, poised_model, 'toplevel');

x = iterate.x;
fx = iterate.feval;
mu = iterate.mu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %%%%%%%%%%%%%%%%%%      END OF DEFT-FUNNEL      %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
