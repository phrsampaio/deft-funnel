function P = deft_funnel_computeP( sampleSet, funY, setting )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Computes the polynomial P, where P is then represented by a row vector
% containing its coefficients for the successive monomials.
% More specifically, these values are:
% P(1)          : constant coefficient,
% P(2:n+1)      : coefficients for the linear terms in x(1)... x(n),
% P(n+2:2n+2)   : coefficients for the squared terms in x(1)^2 ... x(n)^2
% P(2n+3,3n+2)  : coefficients for the quadratic terms of the first subdiagonal:
%                 in x(1)*x(2) ... x(n-1)*x(n)
% P(3n+3,4n+1)  : coefficients for the quadratic terms of the second subdiagonal:
%                 in x(1)*x(3) ... x(n-2)*x(n)
% etc.
%
% Input:
%   - sampleSet : struct of the sample set
%   - funY      : values of the interpolated function at the points in Y
%   - setting   : struct of parameters
%
% Output:
%
%   - P         : a row vector containing the coefficients of the polynomial
%
% Dependencies  : deft_funnel_evalZ
% Called by     : deft_funnel_build_models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ n , p1 ] = size(sampleSet.Y);
q = ( ( n + 1 ) * ( n + 2 ) ) / 2;

if ( setting.whichmodel == 3 && p1 < q )

   % For underdetermined regression model, use min l2-norm model
   setting.whichmodel = 2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( setting.whichmodel == 0 )

   % Build (sub-basis) model (p1 = q) 
   % (sampleSet.QZ and sampleSet.RZ are the factors of Z = M')
   warning off
   P = ( sampleSet.QZ * ( sampleSet.RZ' \ funY' ) )';

elseif ( setting.whichmodel == 1 )

   % Build mixed model: Minimum Frobenius norm model (p1 <= q)
   % and l2-norm model (p1 == n+1 or p1 == q)
   if ( p1 == n+1 || p1 == q )

      warning off
      P(1:p1) = ( sampleSet.RZ \ ( sampleSet.QZ' * funY' ) )';

      if ( p1 == n+1 )
         P(n+2:q) = 0;
      end

   else

      % Build Minimum Frobenius norm model (p1 <= q)
      % minimizing the norm of the Hessian
      % (QZ and RZ are the factors of F = [MQMQ' ML; ML' 0])

      % Compute right-hand side with function values and zero
      rhs = [ funY zeros(1,n+1) ];

      warning off
      mualpha  = ( sampleSet.RZ \ ( sampleSet.QZ' * rhs' ) )';

      % Constant and linear part of P
      P(1:n+1) = mualpha( p1+1: p1+n+1 )';

      % Quadratic part of P
      M        = deft_funnel_evalZ( sampleSet.Y, q )';
      P(n+2:q) = M( :, n+2 : q )' * mualpha( 1: p1 )';

   end

elseif ( setting.whichmodel == 2 )
    
   % Minimum L2 norm model (p1 <= q)
   % (sampleSet.QZ and sampleSet.RZ are the factors of Z = M')
   % Take pseudo-inverse for underdetermined system because the result
   % is different from using backslash-operator
   if ( p1 < q )

      warning off
      P        = ( sampleSet.QZ * ( pinv(sampleSet.RZ') * funY' ) )';

   else

      warning off
      P        = ( sampleSet.QZ * ( sampleSet.RZ' \ funY' ) )';

   end

elseif ( setting.whichmodel == 3 )
    
   % Build Regression model (p1 >= q)
   % (sampleSet.QZ and sampleSet.RZ are the factors of Z = M)
   % Take pseudo-inverse for solving the system because the result
   % is different from using backslash-operator (except for p1==q)
   warning off
   P        = ( pinv(sampleSet.RZ) * sampleSet.QZ' * funY' )';
   
end

end % end of deft_funnel_computeP