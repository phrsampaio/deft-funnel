function [ sampleSet, pos ] = deft_funnel_include_in_Y( sampleSet, x,       ...
    setting, choice_set, poisedness_threshold, criterion, succ )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desc: Attempts to include x in the interpolation set by replacing an existing
% interpolation point. This point is chosen in the set defined by the
% intersection of choice_set and the set of indices such that the absolute
% value of the associated polynomial at x is larger than
% poisedness_threshold. If more than one point remains, the choice is then
% made my maximizing a criterion.  If criterion is 'weighted', then one
% chooses to maximize ||Y(:,j)-x||^2 |L_j(x)|.  If criterion is 'standard',
% then one chooses to maximize |L_j(x)|, while the value of ||Y(:,j)-Y(:,1)||
% alone is used when criterion is 'distance'. (Ties are broken arbitrarily). 
% If no point satisfies the desired properties, nothing is done.
%
% Input: 
%   - sampleSet: structure of the sample set 
%   - x: point to be included into the sample set
%   - setting: parameters' structure
%   - choice_set: a predefined set of indices
%   - poisedness_threshold 
%   - criterion: 'weighted', 'standard' or 'distance' (see desc above)
%   - succ: 1 for successful iteration and 0 otherwise
%
% Output:
%   - sampleSet: update structure of the sample set
%   - pos: index of the replaced interpolation point
%
% Dependencies: deft_funnel_evalL.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y = sampleSet.Y;

% Compute the choice set combining the two requests.
if ( isempty( choice_set ) )
   pos = 0;
   return;
end

% Evaluate the values of the Lagrange polynomials in the choice set at x.
Lvals = deft_funnel_evalL( sampleSet, x, setting, choice_set );

% Compute the combined choice set.
choice = find( abs( Lvals ) > poisedness_threshold );  % the combined choice set
lc     = length( choice );                             % its size

% No suitable point found: return without any action.
if ( lc == 0 )
   pos  =  0;
   return;
end

% Suitable replacement points exist: maximize the criterion.
crit_val = 0;
pos      = 0;
for i = 1:lc
   j = choice( i );
   if ( strcmp(criterion,'weighted') )
      if ( succ == 1 )
         cv = norm( Y(:,j) - x )^2 * abs( Lvals( j ) );
      else
         cv = norm( Y(:,j) - Y(:,1) )^2 * abs( Lvals( j ) );
      end
   elseif ( strcmp(criterion,'standard') )
      cv = abs( Lvals( j ) );
   elseif ( strcmp(criterion,'distance') )
      if ( succ == 1 )
         cv = norm( Y(:,j) - x );
      else
         cv = norm( Y(:,j)-Y(:,1) );
      end
   end
   if ( cv > crit_val )
      pos      = j;
      crit_val = cv;
   end
end

if ( pos == 0 )
    return
end

% Replace Y(:,pos) by x.
sampleSet = deft_funnel_replace_in_Y( sampleSet, x, pos, setting );

end % end of deft_funnel_include_in_Y