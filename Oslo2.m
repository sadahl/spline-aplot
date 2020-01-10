function [bvec] = Oslo2(p, tauvec, cvec, tvec)
%Oslo-Algorithm 2. (Compendium: algorithm 2.20/4.10)
%Knot conversion algorithm:
%Computes B-spline coefficients bvec of f in S_p_t.
%cvec is the original coefficients of f in S_p_tau.
%p is the order of the spline. (p=3: cubic).

m=length(tvec)-p-1;
mu=1;
bvec = zeros(1,m); %preallocation
for i=1:m
    mu=findMu(tauvec, tvec(i), mu);
    if p==0
        bvec(i)=cvec(mu);
    else
        cp=cvec(mu-p:mu);
        offset=p-mu+1; %index offset for cp
        for k = p:-1:1 %algorithm 2.20 modified
            x=tvec(i+k);
            for j = mu:-1:mu-k+1 % overwriting initial cp from the right
                denom = tauvec(j+k)-tauvec(j);
                cp(j+offset) = ( (tauvec(j+k)-x)*cp(j-1+offset) + ...
                                (x-tauvec(j))*cp(j+offset) ) / denom;
            end
        end
        bvec(i) = cp(p+1);
    end % end if
end % end for

end %end function

