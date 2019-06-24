function [ ynew, improvement, msgTR ] = deft_funnel_find_new_yj( sampleSet, ...
             j, setting, Delta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Computes a point which is best to replace yj, the j-th (j>1) column of Y (the
% base for the j-th polynomial) in a ball of radius Delta centered at the
% first column of Y.  This is achieved by maximizing the absolute value of
% the j-th Lagrange polynomial in that ball.  
%
% For conventions on how polynomals are represented, see the documentation of
% evalZ.
%
% Input:
%   - sample       : struct of the new point
%   - j            : position of the sample point in Y to be replaced
%   - setting      : struct of parameters
%   - Delta        : radius
%
% Output:
%
%   - ynew         : the best replacement for Y(:,j)
%   - improvement  : the improvement in poisedness obtained by the update, 
%                    which is equal to |L_j(new y)|. If this value is 
%                    smaller than the threshold input parameter, L and X 
%                    are unchanged by the procedure. 
%   - msgTR        : info message about the Lagrangian polynomials or the
%                    optimization problem involving the Lagragian
%                    polynomials.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Retrieve necessary data
eps_L       = setting.eps_L;
lSolver     = setting.lSolver;

n           = size( sampleSet.Y, 1 );
ynew        = zeros( 1, n );
improvement = 0;

if ( setting.verbose > 2 )
   verbose = 1;
else
   verbose = 0;
end

if ( verbose > 0 )
   disp('--------- enter deft_funnel_find_new_yj ')
end

if ( j < 2 ) % Never attempt to replace the current iterate
   return;
end

% Get the j-th Lagrange polynomial 
Lj = deft_funnel_computeLj( sampleSet, j, setting);

if (~isempty(find(isnan(Lj), 1)) || ~isempty(find( ~isreal(Lj), 1)) || ...
     ~isempty(find(isinf(Lj), 1)))

   msgTR = 'Error0: Lagrange polynomial contains NaN or Inf or nonreal components!!';
   if ( verbose )
      disp(msgTR);
   end
   return
end

% Maximize Lj in a larger 2-norm TR if using infty-norm in the local solver (CG)
if ( lSolver == 2 )
	Delta = sqrt(n)*Delta;
end

% Get the polynomial's gradient and Hessian at the current iterate

% When no shift occurs, the current iterate is Y(:,1)
g  = deft_funnel_gradP( Lj, sampleSet.Y(:,1) );
H  = deft_funnel_hessP( Lj, sampleSet.Y(:,1) );

% Minimize this polynomial and its opposite.
[ pstep, lambda, norms, pvalue, gplus, nfact, neigd, msgTR ] = ...
     deft_funnel_solve_TR_MS( g, H, Delta, eps_L );
[ mstep, lambda, norms, mvalue, gplus, nfact, neigd, msgTR ] = ...
     deft_funnel_solve_TR_MS( -g, -H, Delta, eps_L );

if ( verbose > 0 )
   disp([ ' deft_funnel_find_new_yj: j = ' int2str(j), ' positive value = ', ...
       num2str(pvalue),' step:' ])
   pstep'
   disp([ ' deft_funnel_find_new_yj: j = ' int2str(j), ' negative value = ', ...
       num2str(mvalue),' step:' ])
   mstep'
end

% Select the maximum in absolute value
if ( mvalue < pvalue )
   improvement = abs( mvalue );
   ynew        = sampleSet.Y(:,1) + mstep;
else
   improvement = abs( pvalue );
   ynew        = sampleSet.Y(:,1) + pstep;
end
if ( verbose > 0 )
   disp('--------- exit deft_funnel_find_new_yj ')
end

end % end of deft_funnel_find_new_yj