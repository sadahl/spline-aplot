function [ mu, r, crit ] = C_lengthratio( cpoints, lsize, ~ )
% [ mu, r, crit ] = C_lengthratio cpoints,lsize,~  )
%   Determines the index mu of the control point assumed to hold greatest
%   visual improvement potential.
%   lsize sets the nb of points in local control polygon (lcp). Default: 3
%       This parameter remains to study closer.
%   The range parameter r=0 is kept for compatibility with other methods.
%   The criterion crit is the greatest ratio between the length of the
%   local control polygon and the length of the local base line.
%   The index mu of an (the) inner control point is kept.
%   If more than one point satisfies the criterion, 
%   the first (counted from the left) is returned.

r=0;
n=length(cpoints);

if nargin<=2
    lsize = 3; % one inner point per local control polygon
end
crit=0;

% Using simplified length-expressions
% diffs = abs( cpoints(:,2:n)-cpoints(:,1:n-1) );
% lengths = max(diffs) + 1/2*min(diffs); % the approximation of norm (s=1,2)
%     % lengths(i) is now approx length of [cp_{i} cp_{i+1}]
% idxs = 1:n-lsize+1; % indices of start and end of each lcp
% idxe = lsize:n;
% baselinediffs = abs( cpoints(:,idxe)-cpoints(:,idxs) );
% baselinelengths = max(baselinediffs) + 1/2*min(baselinediffs);

% Alternative, using exact length
diffs = cpoints(:,2:n)-cpoints(:,1:n-1);
lengths = sqrt( dot(diffs, diffs) );
%     % lengths(i) is now approx length of [cp_{i} cp_{i+1}]
idxs = 1:n-lsize+1; % indices of start and end of each lcp
idxe = lsize:n;
baselinediffs = cpoints(:,idxe)-cpoints(:,idxs);
baselinelengths = sqrt( dot(baselinediffs, baselinediffs) );

% We compute the partial sums of lengths corresponding to each lcp
% Improve this implementation, for instance using sum of matrix rows
for i=1:n-lsize+1
%     lengths(i:i+lsize-2)
    lcplengths(i)=sum (lengths(i:i+lsize-2) );
end
% lcplengths
candidates = lcplengths./baselinelengths; % The length-ratios
crit = max(candidates);
mu = find(~(candidates-crit),1) + 1; % adjusts for right-point association
                                     % makes sense only for lsize=3
                                     % consider adjusting using r later.

end


%% Testing
% c =
% 
%     6.8505    6.7839    6.7735    2.7460    1.6455
%     8.6326    7.0633    0.4289    8.1367    2.4993
% 
% C_lengthratio(c)
% 
% baselinelengths =
% 
%     8.2040    4.1782    5.5302
% 
% 
% lcplengths =
% 
%     8.2051   15.3310   14.4405
% 
% 
% candidates =
% 
%     1.0001    3.6693    2.6112
% 
% 
% ans =
% 
%      3