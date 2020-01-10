function xm = midpointRefine(x)
%midpointRefine gives the refined vector/grid of vector/grid x by inserting
%midpoints. Returns only the new points
%   x will typically be a vector/grid of knot averages and thus will not
%   contain duplicates.
%   1D: function returns a 1D vector containing one
%   point in the middle of each interval. Input assumed to be a sorted
%   row-vector. Output is also sorted.

    xuniq=knt2brk(x); %the unique knots
    xm = (xuniq(1:end-1)+xuniq(2:end))/2; %all the midpoints
    xm=sort([x xm]);

end

