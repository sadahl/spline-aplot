function [ mu, r, crit ] = C_baselineDist( cpoints, lsize, ~ )
% [ mu, r, crit ] = C_baselineDist cpoints,~,~  )
%   Determines the index mu of the control point assumed to hold greatest
%   visual improvement potential.
%   lsize sets the nb of points in local control polygon. Default: 3
%       This parameter remains to study closer.
%   The range parameter r=0 is kept for compatibility with other methods.
%   The criterion crit is the greatest distance to local control polygon
%   base line.
%   If more than one point satisfies the criterion, 
%   the first (counted from the left) is returned.
%   Implementation can be improved by vectorizing the called subroutine
%   baselineDisti.

r=0;
n=length(cpoints);

if nargin<=2
    lsize = 3; % one inner point per local control polygon
end
crit=0;

for i=1:n-lsize+1 % each local control polygon
    for j=i+1:i+lsize-2 % each inner point in the current local cpolygon
        candidate = baselineDisti(  cpoints(:,i), ...
                                    cpoints(:,i+lsize-1), ...
                                    cpoints(:,j)  );
        if candidate > crit
            crit = candidate;
            mu = j;
        end
    end
end

end


%% Testing
% cc =
% 
%     8.5434    9.8915    3.3346    3.6097    2.5168
%     0.1951    5.6788    6.4881    6.0683    0.2846
% 
% C_baselineDist(cc)
% 
% i =
% 
%      1
% 
% 
% j =
% 
%      2
% 
% 
% candidate =
% 
%    20.5670
% 
% 
% i =
% 
%      2
% 
% 
% j =
% 
%      3
% 
% 
% candidate =
% 
%     0.2519
% 
% 
% i =
% 
%      3
% 
% 
% j =
% 
%      4
% 
% 
% candidate =
% 
%     0.1073
% 
% 
% ans =
% 
%      2


