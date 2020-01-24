function [ bool ] = sameSide( cs,ce,pi,pj )
%[ bool ] = sameSide( cs,ce,pi,pj )
%   Checks whether points pi, pj are on the same side of line (cs,ce).
%   Returns true if pi and pj are on the same side
%   Returns false if not or if one of the points is exactly on  (cs,ce).

dpi=pi-cs;
dpj=pj-cs;
dp=ce-cs;

% cross requires dimension 3
if length(pi)==2
    dpi(:,3)=0;
    dpj(:,3)=0;
    dp(:,3)=0;
end

bool = ( dot(cross(dpi,dp),cross(dpj,dp))>0 );

end

