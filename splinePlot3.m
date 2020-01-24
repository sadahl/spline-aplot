function spr =  splinePlot3(sp, spr, newt)
%Plots the spline sp and the control polygons associated with the spline
%representations.
% The old spline is sp, the new spline is spr (refined).
% The knots that were inserted are held in array newt.

% Setting up the plot window
p=sp.order-1; % order of spline
scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/1.3 scrsz(3)/1.3 scrsz(4)/1.3])

% Plotting control polygon,
% according to the dimensions of the spline coefficients
if sp.dim==1
    c1=aveknt(sp.knots,p+1);
    c2=sp.coefs;
    c1r=aveknt(spr.knots,p+1);
    c2r=spr.coefs;    
else % sp.dim==2
    c1=sp.coefs(1,:);
    c2=sp.coefs(2,:);
    c1r=spr.coefs(1,:);
    c2r=spr.coefs(2,:);
end % end if

pc = plot(c1, c2,'--ok');%, 'DisplayName', legcp);  % Old in black
hold on;
pcr = plot(c1r, c2r, '-b'); % New in blue.

if sp.dim==1 && nargin==3
    knotPlotter(sp.knots, sp.coefs, newt);    
end

fnpltpoints=fnplt(sp); %plotted with reference tool fnplt
ps=plot(fnpltpoints(1,:),fnpltpoints(2,:), 'r-'); hold off;
% N=length(sp.coefs);

end


