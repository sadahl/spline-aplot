function [ mu, r, crit ] = C_angledot_simple( cpoints, ~, ~ )
% [ mu, r, crit ] = C_angledot_simple( cpoints,~,~  )
%   Determines the index mu of the control point assumed to hold greatest
%   visual improvement potential.
%   The range parameter r=0 is kept for compatibility with other methods.
%   The criterion crit is the biggest angle of change.
%   Uses an approximation of the angle of directional change
%   from c_{\mu}-c_{\mu-1} to c_{\mu+1}-c_{\mu}.
%   This angle is associated with c_{\mu}.
%   Unit: radians.
%   Uses the transformation T(x) = -pi/2*x*abs(x)+pi/2 to mimic arccos.
%   If more than one point satisfies the criterion, 
%   the first (counted from the left) is returned.

r=0;
n=length(cpoints);

dp = cpoints(:,2:n) - cpoints(:,1:n-1);
dpnormsq = dot( dp(:,1:n-1), dp(:,1:n-1) );

% dot product formula, see thesis section 3.2.2
argnum = dot( dp(:,1:n-2), dp(:,2:n-1) );
argdenomsq = dpnormsq(1:n-2) .* dpnormsq(2:n-1);
% angles = acos(argnum/sqrt(argdenomsq)); % the exact angle formula
% Applying the approximation instead - the transformation T:
angles = -pi/2*argnum.*abs(argnum)./argdenomsq + pi/2;

% Must handle case of subsequent control points being equal:
% If points with index i and i+1 are equal,
% then the angles at these points are undefined.
% The corresponding value with index i in dp is 0.
ind = find(argnum==0);    
angles(ind) = 0; % setting 'undefined' angles to 0 degrees

crit = max(angles);
% NB: this angle value must be compared to a modified angle threshold,
% see ChapterX attached to thesis.
mu = find(~(angles-crit),1) + 1; % adjusts for right-point association

end


%% Testing
% c =
% 
%     5.2019    9.2223    8.4992    9.0958    8.1680
%     0.9591    0.2217    0.5818    4.4719    4.4272
% 
% C_angledot_simple(c)
% 
% dp =
% 
%     4.0204   -0.7231    0.5966   -0.9278
%    -0.7374    0.3601    3.8901   -0.0447
% 
% 
% dpnormsq =
% 
%    16.7072    0.6525   15.4885    0.8628
% 
% 
% angles =
% 
%     3.0211    1.4247    1.6330
% 
% 
% ans =
% 
%      2

