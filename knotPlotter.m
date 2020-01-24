function h = knotPlotter(t, c, newt)
%Plotting knots and their multiplicities. (functional case)
%Knots are shown as red circles slightly beneath the spline.
%t are the knots, shown in circles
%c are the coefficients of the spline
%newt, if provided, are the knots that were inserted, shown as diamonds

ymin=min(c); ymax=max(c); yamp=ymax-ymin;
ybaseline=ymin-yamp*.05; % underneath the lower part of cp
m = knt2mlt(t);
y=ybaseline-m*yamp*.025;
h = scatter(t, y, 'ro', 'filled', 'DisplayName', 'Knot placement');
hold on;

if nargin==3
    newt = sort(newt);
    newm = knt2mlt(newt);
    %shifting vertically when there are multiplicities
    tshare=sort(intersect(newt,t));
    for i=1:length(tshare)
        element=tshare(i);
        nbinnewt=sum(~(newt-element));
        nbint=sum(~(t-element));
        myindex=find(newt==element,1); 
        % first index in newm corresponding to tdiff(i)
        for j=myindex:myindex+nbinnewt-1
            newm(j)=newm(j)+nbint;
        end
    end

    newy=ybaseline-newm*yamp*.025;
    h = scatter(newt, newy, 'md', 'filled', 'DisplayName', 'Knot placement');
end %end if

end

