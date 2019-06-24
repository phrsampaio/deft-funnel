function Lj = deft_funnel_computeLj( sampleSet, j, setting )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Computes the coefficients of the j-th Lagrange polynomial.
%
% Input:
%   - sampleSet : the structure of the sample set
%   - j         : the index of the Lagrange polynomial to be computed
%   - setting   : the structure of the set of parameters
%
% Output:
%
%   - Lj        : a row vector containing the coefficients of the polynomial
%
% Dependencies  : deft_funnel_evalZ
% Called by     : deft_funnel_find_new_yj, deft_funnel_find_new_yj_bc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y = sampleSet.Y;
[ n , p1 ] = size(Y);

% Define the size of the polynomial
if ( setting.whichmodel == 0 )
   q = p1;
else
   q = ( ( n + 1 ) * ( n + 2 ) ) / 2;
end

if ( setting.whichmodel == 3 && p1 < q )

   % For underdetermined regression model, use min l2-norm model
   setting.whichmodel = 2;

end

if ( setting.whichmodel == 0  )

   % Build (sub-basis) Lagrange polynomials (p1 = q)
   % (QZ and RZ are the factors of Z = M')
   warning off
   Lj = ( sampleSet.QZ * ( sampleSet.RZ' \ [zeros(j-1,1);1;zeros(p1-j,1)]) )';

elseif ( setting.whichmodel == 1 )
    
   % Build mixed Lagrange polynomials:
   % Minimum Frobenius norm Lagrange polynomials (p1 <= q) or
   % L2-norm Lagrange polynomials (p1 == n+1 or p1 == q )
   if ( p1 == n+1 || p1 == q )

      % Build linear/quadratic Lagrange polynomials
      warning off
      Lj = ( sampleSet.RZ \  sampleSet.QZ' * [zeros(j-1,1);1;zeros(p1-j,1)] )';

      if ( p1 == n+1 )
         Lj(n+2:q) = 0;
      end

   else

      % Build Minimum Frobenius-norm Lagrange polynomials (p1 <= q)

      % Compute the right-hand side of the system
      rhs = [zeros(j-1,1);1;zeros(p1+n+1-j,1)] ;

      warning off
      mualpha    = ( sampleSet.RZ \ ( sampleSet.QZ' * rhs ) )';

      % Constant and linear part of P
      Lj(1:n+1) = mualpha( p1+1: p1+n+1 )';

      % Quadratic part of P
      M         = deft_funnel_evalZ( sampleSet.Y, q )';
      Lj(n+2:q) = M( :, n+2 : q )' * mualpha( 1: p1 )';

   end

elseif ( setting.whichmodel == 2 )
    
   % Build Minimum L2 norm Lagrange polynomials (p1 <= q)
   % (QZ and RZ are the factors of Z = M')
   % Take pseudo-inverse for underdetermined system because the result
   % is different from using backslash-operator
   if ( p1 < q )

      warning off
      Lj = ( sampleSet.QZ * ( pinv(sampleSet.RZ') * [zeros(j-1,1);1;zeros(p1-j,1)]) )';

   else

      warning off
      Lj = ( sampleSet.QZ * ( sampleSet.RZ' \ [zeros(j-1,1);1;zeros(p1-j,1)]) )';

   end

elseif ( setting.whichmodel == 3 )
    
   % Build Regression Lagrange polynomials (p1 >= q)
   % (QZ and RZ are the factors of Z = M)
   warning off
   Lj = ( pinv(sampleSet.RZ) *  sampleSet.QZ' * [zeros(j-1,1);1;zeros(p1-j,1)] )';

end

end % end of deft_funnel_computeLj