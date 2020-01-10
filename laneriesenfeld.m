function [c1n,c2n] = laneriesenfeld(c1,c2,first)
%Implementens one round of the subdivision method of Lane and Riesenfeld
%   Takes as input
%   a 2D control polygon (c1,c2) (c1 may be the knot averages tstar)
%   c1 and c2 are assumed to be of same length, n.
%   Optional argument first to indicate only to run the first round
%
%   Inserts midpoints in first round.
%   Cuts the corners of the control polygon at half the length of
%   adjoining segments.
%   Returns the new control polygon, which better approximates a smooth
%   curve.

n=length(c1);

if nargin==3 && first
    %first round
    %allocating
    c1n=zeros(1,2*n-1);
    c2n=zeros(1,2*n-1);
    %computing new values
    for i=1:n-1
        c1n(2*i-1)=c1(i);
        c2n(2*i-1)=c2(i);
        
        c1n(2*i)=( c1(i) + c1(i+1) )/2;
        c2n(2*i)=( c2(i) + c2(i+1) )/2;
    end %end for
    %last index
        c1n(2*n-1)=c1(n);
        c2n(2*n-1)=c2(n);
    
else
    %general round
    %allocating
    c1n=zeros(1,n-1);
    c2n=zeros(1,n-1);
    %computing new values
    for i = 1:n-1
        c1n(i)=( c1(i) + c1(i+1) )/2;
        c2n(i)=( c2(i) + c2(i+1) )/2;
    end %end for
end %end if
end

