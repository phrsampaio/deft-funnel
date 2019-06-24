function [ sampleSet, iterate ] = deft_funnel_set_current_subspace( ...
    oldSampleSet, oldIterate, vstatus, const)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Updates iterate and sampleSet structures for the current subspace.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iterate   = oldIterate;
sampleSet = oldSampleSet;

iterate.indfree   = find( vstatus == const.free );
iterate.indfix    = find( vstatus >= const.fixed );
iterate.x         = sampleSet.X( iterate.indfree, sampleSet.i_xbest );
iterate.nfix      = length( iterate.indfix );
iterate.xdim      = length( iterate.indfree );
iterate.feval     = sampleSet.fX( sampleSet.i_xbest );
iterate.ceval     = sampleSet.cX( :, sampleSet.i_xbest );

sampleSet.Y       = sampleSet.X( iterate.indfree, sampleSet.ind_Y );
sampleSet.fY      = sampleSet.fX( sampleSet.ind_Y );
sampleSet.cY      = sampleSet.cX( :, sampleSet.ind_Y );

end % enf of deft_funnel_set_current_subspace