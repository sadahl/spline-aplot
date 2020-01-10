function [ xx ] = knotInteriors( sp )
%Returns the knot interiors of the spline sp
% Part of insertFunction family used by splinePlot2

xx=setdiff(sp.knots,[sp.knots(1) sp.knots(sp.number+sp.order)])
% xx= [xx xx xx]

end

