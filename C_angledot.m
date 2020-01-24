function [ mu, r, crit ] = C_angledot( cpoints, ~, ~ )
% [ mu, r, crit ] = C_angledot( cpoints,~,~  )
%   Determines the index mu of the control point assumed to hold greatest
%   visual improvement potential.
%   The range parameter r=0 is kept for compatibility with other methods.
%   The criterion crit is the biggest angle of change.
%   Computes the angle of directional change
%   from c_{\mu}-c_{\mu-1} to c_{\mu+1}-c_{\mu}.
%   The greatest such angle is associated with c_{\mu}.
%   Unit: radians.
%   If more than one point satisfies the criterion, 
%   the first (counted from the left) is returned.

r=0;
n=length(cpoints);

dp = cpoints(:,2:n) - cpoints(:,1:n-1)
dpnormsq = dot( dp(:,1:n-1), dp(:,1:n-1) )

% dot product formula, see thesis section 3.2.2
argnum = dot( dp(:,1:n-2), dp(:,2:n-1) );
argdenomsq = dpnormsq(1:n-2) .* dpnormsq(2:n-1);

% Must handle case of subsequent control points being equal:
% If points with index i and i+1 are equal,
% then the angles at these points are undefined.
% The corresponding value with index i in dp is 0.
ind = find(argnum==0);    
% arg = argnum./sqrt(argdenomsq);  % risky for floating point calc.
arg = sqrt(argnum.^2./argdenomsq) % may hold division by zero, here handled
arg(ind) = 1; % setting 'undefined' angles to 0 degrees

% Collecting other cases, as a floating-point failsafe
ind = [find(arg<-1) find(arg>1)];
if ~isempty(ind)
    disp('Warning: angle calculation truncated, invalid arg for arccos')
    arg(ind) = -1;
end
angles = acos(arg); % the exact angle formula

crit = max(angles);
mu = find(~(angles-crit),1) + 1; % adjusts for right-point association

end


%% Testing
% c =
% 
%     3.1332    7.5259    4.0165    0.2338    5.3283
%     1.4645    3.2527    0.1438    8.7576    2.4761
% 
% C_angledot(c)
% 
% angles =
% 
%     2.8032    1.8820    2.8740
% 
% 
% ans =
% 
%      4