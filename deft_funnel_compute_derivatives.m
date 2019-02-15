function derivatives = deft_funnel_compute_derivatives( models, sampleSet, iterate )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Desc: Computes the derivatives of the surrogate models.
%
% Input:
%   - models      : struct of the surrogate models
%   - Ynew        : struct of the sample set
%   - setting     : struct of pthe current iterate
%
% Output:
%   - derivatives : struct of the derivatives at the current iterate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gradient and hessian of the objective function
gfx = deft_funnel_gradP( models.f, iterate.x );
Hf  = deft_funnel_hessP( models.f, iterate.x );

% Jacobian of the contraint functions
m = size( sampleSet.cY, 1 );
for i=1:m
    J(i,:) = deft_funnel_gradP( models.c(i,:), iterate.x );
end

% Hessian matrices of the constraint functions
HClist = [];
for i=1:m
    Ci = deft_funnel_hessP( models.c(i,:), iterate.x );
    HClist = [ HClist Ci ];
end

% Assemble them all into a single structure
derivatives.gfx = gfx;
derivatives.Hf = Hf;
derivatives.J = J; 
derivatives.HClist = HClist;

end % enf of deft_funnel_compute_derivatives
