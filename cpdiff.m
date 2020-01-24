function [ cdiff ] = cpdiff( cpoints,cpointsnew,mu,p )
% [ cdiff ] = cpdiff( cpoints,cpointsnew,mu,p )
%   Computes a measure of the difference of the new and old control polygon
%   Computes the maximum change in coordinates of the control points which
%   are 'moved'.
%   Possible improvement: extract this from the knot insertion function.
%
%   cpoints are the old control points
%   cpointsnew are the new control points
%   mu is the index such that the inserted knot z lies in 
%   the interval [t_mu^*, t_{mu+1}^*).
%   p is the degree of the spline

cdiff=0;

if nargin == 4
    % Method 1: use input information and Bohms equations to locate the
    % coordinates that have changed
    for i=mu-p+1:mu
        cdiff=max([norm(cpoints(:,i-1)-cpointsnew(:,i));...
                   norm(cpoints(:,i)  -cpointsnew(:,i));
                   cdiff]);
        %cheaper alternative:
        %cdiff=max([length_simple(cpoints(:,i-1),cpointsnew(:,i));
        %           length_simple(cpoints(:,i)  ,cpointsnew(:,i));
        %           cdiff])
               
    end
else
    % Method 2: search for changes, from left to right
    % in the coordinate vectors
    i=1;
    tol=eps(1);
    while max(abs(cpoints(:,i)-cpointsnew(:,i))) < tol
        i=i+1;
    end
    i=i-1;
    while max(abs(cpoints(:,i)-cpointsnew(:,i+1))) > tol
        i=i+1;
        cdiff=max([norm(cpoints(:,i-1)-cpointsnew(:,i));...
                   norm(cpoints(:,i)  -cpointsnew(:,i));
                   cdiff]);

    end
end %end if

