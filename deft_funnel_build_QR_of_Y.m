function [ sampleSet, badcond ] = deft_funnel_build_QR_of_Y( sampleSet, setting )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Desc: Computes the QR factorization of the matrix containing
%  the polynomial expansion of the interpolation points. 
%
%  Input:
%    - sampleSet : struct of the sample set
%    - setting   : struct of parameters
%
%  Output:
%    - sampleSet : struct of the sample set containing the QR factors of the
%                  matrix containing the polynomial expansion of the 
%                  interpolation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y = sampleSet.Y;
[ n , p1 ] = size(Y);

% Define the tolerance for ill conditioning
tol = 1.0e-10;
badcond = 0;

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

if ( setting.whichmodel == 0 )

   % QR of (sub-basis) interpolation matrix (p1 = q) or
   % (factorize matrix Z = M')
   Z = deft_funnel_evalZ( Y, q );

   % Check condition of Z and cure if ill-conditioned
   if ( isempty(find(isnan(Z), 1)) && isempty(find(isinf(Z), 1)) )
      condZ = cond(Z);
      if (condZ > setting.kappa_ill)
         badcond = 1;
      end
   else
      badcond = 1;
   end

   if ( badcond )
      [U,S,V]        = svd(Z,0);
      Sdiag          = diag(S);
      indices        = find(Sdiag < tol);
      Sdiag(indices) = tol;
      %S(eye(size(S))~=0) = Sdiag;
      S              = diag(Sdiag);
      
      M              = (V * S * U')';
      [ QZ, RZ ]     = qr( M );
   else
      [ QZ, RZ ]     = qr( Z );
   end
   
elseif ( setting.whichmodel == 1 )

   % Mixed model: minimum Frobenius-norm model (when underdetermined) and 
   % minimum l2-norm model (at linear and quadratic degree)
   if ( p1 == n+1 || p1 == q )

      F = deft_funnel_evalZ( Y, p1 )';

   else

      % QR of Minimum Frobenius norm interpolation matrix (p1 <= q)
      % (factorize matrix F = [MQMQ' ML; ML' 0])
      M  = deft_funnel_evalZ( Y, q )';
      ML = M( :, 1:n+1 );
      MQ = M( :, n+2:q );
      F  = [ MQ * MQ' ML; ML' zeros(n+1,n+1) ];

   end

   %  Check condition of Z and cure if ill-conditioned
   if ( isempty(find(isnan(F), 1)) && isempty(find(isinf(F), 1)) )
      condZ = cond(F);
      if (condZ > setting.kappa_ill)
         badcond = 1;
      end
   else
      badcond = 1;
   end

   if ( badcond )
      [U,S,V]        = svd(F,0);
      Sdiag          = diag(S);
      indices        = find(Sdiag < tol);
      Sdiag(indices) = tol;
      S              = diag(Sdiag);
      M              = (V * S' * U')';
      [ QZ, RZ ]     = qr( M );
   else
      [ QZ, RZ ] = qr( F );
   end

elseif ( setting.whichmodel == 2 )

   % QR of Minimum L2 norm interpolation matrix (p1 <= q)
   % (factorize matrix Z = M')
   Z = deft_funnel_evalZ( Y, q );

   % Check condition of Z and cure if ill-conditioned
   if ( isempty(find(isnan(Z), 1)) && isempty(find(isinf(Z), 1)) )
      condZ = cond(Z);
      if (condZ > setting.kappa_ill)
         badcond = 1;
      end
   else
      badcond = 1;
   end

   if ( badcond )
      [U,S,V]        = svd(Z,'econ');
      Sdiag          = diag(S);
      indices        = find(Sdiag < tol);
      Sdiag(indices) = tol;
      S              = diag(Sdiag);
      M              = (V * S' * U')';
      [ QZ, RZ ]     = qr( M );
   else
      [ QZ, RZ ]     = qr( Z );
   end

elseif ( setting.whichmodel == 3 )

   % QR of Regression interpolation matrix (p1 >= q)
   % (factorize matrix Z = M)
   Z =  deft_funnel_evalZ( Y, q )' ;

   % Check condition of Z and cure if ill-conditioned
   if ( isempty(find(isnan(Z), 1)) && isempty(find(isinf(Z), 1)) )
      condZ = cond(Z);
      if (condZ > setting.kappa_ill)
         badcond = 1;
      end
   else
      badcond = 1;
   end

   if ( badcond )
      [U,S,V]        = svd(Z,0);
      Sdiag          = diag(S);
      indices        = find(Sdiag < tol);
      Sdiag(indices) = tol;
      S              = diag(Sdiag);
      M              = (V * S * U')';
      [ QZ, RZ ]     = qr( M );
   else
      [ QZ, RZ ]     = qr( Z );
   end

end

sampleSet.QZ = QZ;
sampleSet.RZ = RZ;

end % end of deft_funnel_build_QR_of_Y
