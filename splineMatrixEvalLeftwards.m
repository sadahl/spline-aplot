function [fx, mu] = splineMatrixEvalLeftwards(p, tvec, c0, x, muin)
%Obsolete: use matlab native function fnval(sp,x) where sp is the spline.
%
%Computes the function value of a spline.
%p is the order of the spline (p=3: cubic)
%tvec is the knot vector
%c0 is the vector holding the B-spline coefficients
%x is the evaluation point
%muin is an optional argument that may speed up the computation

if nargin == 5
    mu=findMu(tvec, x, muin);
else
    mu=findMu(tvec, x);
end

for k = p:-1:1
    for j = mu:-1:mu-k+1 % overwriting initial c0 from the right
        denom = tvec(j+k)-tvec(j);
        c0(j) = ( (tvec(j+k)-x)*c0(j-1) + (x-tvec(j))*c0(j) ) / denom;
    end
end
fx = c0(mu);
end % end function

