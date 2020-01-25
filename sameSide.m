function [ bool ] = sameSide( cs,ce,p1,p2 )
%[ bool ] = sameSide( cs,ce,pi,pj )
%   Checks whether points pi, pj are on the same side of line (cs,ce).
%   Returns true if pi and pj are on the same side
%   Returns false if not or if one of the points is exactly on  (cs,ce).

dp1=p1-cs;
dp2=p2-cs;
dp=ce-cs;

% cross requires dimension 3
if length(p1)==2
    dp1(:,3)=0;
    dp2(:,3)=0;
    dp(:,3)=0;
end

bool = ( dot(cross(dp1,dp),cross(dp2,dp))>0 );

end

