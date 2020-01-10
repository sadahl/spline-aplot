function [ dSq ] = baselineDist(cs,ce,p)
%Returns the square of the distance d from a point to a line segment.
%   Uses dot products to compute distances.
%   If the projection is outside the segment between the points,
%   the method returns the distance to the closest point.
%   cs is the start point of the line segment
%   ce is the end point of the line segment
%   p is the specified point to which we measure distance.

%Moving the origin to cs
x2=ce(1)-cs(1);
y2=ce(2)-cs(2);
px=p(1)-cs(1);
py=p(2)-cs(2);

dotprod = px*x2 + py*y2;
% projlenSq;
if dotprod <= 0.0
    %the projection of p is outside the segment, on side of cs
    projlenSq=0;
else
    %moving origin to ce for similar check
    px=x2-px;
    py=y2-py;
    dotprod=px*x2 + py*y2;  % a dot b = -a dot -b
    if dotprod <= 0.0
        projlenSq = 0.0;
    else
        %p is between cs and ce
        projlenSq = dotprod * dotprod / (x2*x2 + y2*y2);
    end
end
dSq = max(0,px*px + py*py - projlenSq);

end