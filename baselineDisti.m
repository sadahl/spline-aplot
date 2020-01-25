function [ dSq ] = baselineDisti(cs,ce,p)
%[ dSq ] = baselineDist(cs,ce,p)
% Returns the square of the distance d from a point to a line segment.
%   Improved version of baselineDist, using MATLAB built-in functions.
%   Uses dot products to compute distances.
%   If the projection is outside the segment between the points,
%   the method returns the distance to the closest point.
%   cs is the start point of the line segment
%   ce is the end point of the line segment
%   p is the specified point to which we measure distance.

%Moving the origin to cs
bs=ce-cs;
pt=p-cs;

dotprod = dot(bs,pt);
% projlenSq;
if dotprod <= 0.0
    %the projection of p is outside the segment, on side of cs
    projlenSq=0;
else
    %moving origin to ce for similar check
    pt=bs-pt;
    dotprod=dot(bs,pt);  % a dot b = -a dot -b
    if dotprod <= 0.0
        projlenSq = 0.0;
    else
        %p is between cs and ce
        projlenSq = dotprod * dotprod / norm(bs)^2;
    end
end
dSq = max(0,norm(pt)^2 - projlenSq);

end