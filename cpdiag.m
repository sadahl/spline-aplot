function [ cdelta ] = cpdiag( c1,c2 )
%Computes the diagonal of a control polygon
%   Uses max/min-tests to find the encompassing rectangle containing the
%   control polygon. Returns the diagonal length of this rectangle.
ll=[min(c1),min(c2)];
ur=[max(c1),max(c2)];
cdelta=norm(ur-ll);

end

