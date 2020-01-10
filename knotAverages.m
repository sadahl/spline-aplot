function xx = knotAverages(sp)
%Returns the knot averages of the spline sp, excluding the end points
% Part of insertFunction family used by splinePlot2

xx = aveknt(sp.knots,sp.order);
xx = setdiff(xx,[sp.knots(1) sp.knots(sp.number+sp.order)]);

end

