function H = deft_funnel_buildHLag( mu, derivatives )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Builds the hessian of the Lagrangian function.
% 
% Input:
%   - mu          : Lagrange multipliers
%   - derivatives : struct of derivatives
%
% Output: 
%   - H           : Hessian matrix of the Lagrangian function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = size( mu, 1 );
H = derivatives.Hf;
n = size( derivatives.HClist, 1 );

for i = 1:m
    H = H + mu( i )*derivatives.HClist( 1:n, (i-1)*n+1:i*n );
end

end % end of deft_funnel_buildHLag