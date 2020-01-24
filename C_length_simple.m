function [ mu, r, crit ] = C_length_simple( cpoints, ~, ~ )
% [ mu, r, crit ] = C_length_simple( cpoints,~,~  )
%   Determines the index mu of the control point assumed to hold greatest
%   visual improvement potential.
%   The range parameter r=0 is kept for compatibility with other methods.
%   Segment [c_{\mu-1},c_{\mu}) is associated with c_{\mu}.
%   The criterion crit is an approximation of the longest
%   control polygon segment. If more than one segment satisfies the
%   criterion, the first (counted from the left) is returned.

r=0;
n=length(cpoints);
diffs = abs( cpoints(:,2:n)-cpoints(:,1:n-1) );
lengths = max(diffs) + 1/2*min(diffs); % the approximation of norm (s=1,2)
crit = max(lengths);
mu = find(~(lengths-crit),1) + 1; % adjusts for right-point association

end