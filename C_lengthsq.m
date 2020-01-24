function [ mu, r, crit ] = C_lengthsq( cpoints, ~, ~ )
% [ mu, r, crit ] = C_lengthsq( cpoints,~,~  )
%   Determines the index mu of the control point assumed to hold greatest
%   visual improvement potential.
%   The range parameter r=0 is kept for compatibility with other methods.
%   Segment [c_{\mu-1},c_{\mu}) is associated with c_{\mu}.
%   The criterion crit is the square length of the longest
%   control polygon segment. If more than one segment satisfies the
%   criterion, the first (counted from the left) is returned.

r=0;
n=length(cpoints);
diffs = cpoints(:,2:n)-cpoints(:,1:n-1);

lengthssq = dot(diffs, diffs);
crit = max(lengthssq);
mu = find(~(lengthssq-crit),1) + 1; % adjusts for right-point association

end

