function [c1n,c2n] = chaikin(c1,c2)
%Implements one round of Chaikins algorithm.
%   Takes as input
%   a 2D control polygon (c1,c2) (c1 may be the knot averages tstar)
%   c1 and c2 are assumed to be of same length, n.
%
%   Cuts the corners of the control polygon at one quarter of the length of
%   adjoining segments.
%   Returns the new control polygon, which better approximates a smooth
%   curve.

%allocating
n=length(c1);
c1n=zeros(1,2*n-2);
c2n=zeros(1,2*n-2);

%computing new values
for i = 1:n-1
    c1n(2*i-1)=(3*c1(i)+c1(i+1))/4;
    c2n(2*i-1)=(3*c2(i)+c2(i+1))/4;

    c1n(2*i)=(c1(i)+3*c1(i+1))/4;
    c2n(2*i)=(c2(i)+3*c2(i+1))/4;

end %end for

end

