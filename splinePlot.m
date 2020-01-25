function [spr, splineValuesnew, N] = ...
    splinePlot(sp, maxN, r, method, mytitle, showSpline, showSampledSpline)
%DEPRECATED - KEPT FOR ARCHIVE and
%CONTINUED SUPPORT FOR OLD SCRIPTS ('introPlotting.m')    
%Plots the spline sp and its control polygon. The control polygon is
%refined in 'r' rounds using method 'method'.
%Valid methods: uniform, newknt, optknt, chbknt.
%
%Each round is calculated and displayed with a mouse click.
%The maximum number of plot points allowed is maxN.
%The spline sp must be in B-form.
%The spline itself is plotted when showSpline is True.
%
%The function returns the spline spr in the finer space, the sample values
%splineValuesnew of the refined spline taken at the new knot averages.
%A plot of these points is given if showSampledSpline is
%True.

% Setting up the plot window
scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/1.3 scrsz(3)/1.3 scrsz(4)/1.3])
bufrate=.1; %for second axis adjustment
k=sp.order; p=k-1;
ymin=min(sp.coefs); ymax=max(sp.coefs);
amp=max(2/bufrate,ymax-ymin);
% axis([min(sp.knots) max(sp.knots) ymin-bufrate*amp ymax+bufrate*amp]);
% title({mytitle,method}); %two-liner in title
title(mytitle);
hold on;

% Option for rendering the spline using Matlab built-in method
% pref = plot(xx,fnval(sp,xx),'-','LineWidth',1,'DisplayName','The Spline'); 
if showSpline
    fnpltpoints=fnplt(sp); %plotted with reference tool fnplt
    legsp='The spline';
    ps=plot(fnpltpoints(1,:),fnpltpoints(2,:), 'r.');%,'DisplayName','The Spline');
end

% The Control Polygon:
xxstar=aveknt(sp.knots,p+1); % initial sampling at the knot averages
N=length(xxstar);
legcp='Control Polygon';
pc = plot(xxstar, sp.coefs,'-ok');%, 'DisplayName', legcp); 
% legend('','Position', [.5 .5 .5 .5]);
% Comparison: sampling the spline at the knot average locations
if nargin==7 && showSampledSpline
    mystr=sprintf('Sampling, N=%d',N);
    splineValues=fnval(sp,xxstar);
    psample=plot(xxstar,splineValues,'m-.o','DisplayName',mystr);
end

% Preparing refinement rounds
mylegend=sprintf('Refined control Polygon, N=%d',N);
pr = plot(xxstar, sp.coefs,'-');%, 'DisplayName', mylegend); 
% legend('show');
   
% Choosing how to refine the knot vector/the control polygon
origxx=sp.knots; % using the original knots, not the knot averages
xx=origxx;
switch method
    case 'uniform' %Uniform distribution within each knot interval
         % Doubling the number of points (adding midpoints)
         % We may plot too few points (doubling is rough)
         % Note: see built-in function syntax: sprfn(sp)
        while N<=maxN/2 %not precise: what is the number of knot averages?
            xx=midpointRefine(xx);
            N=length(xx);
        end
        
    case 'newknt' %Old Fortran method, see matlab documentation
        xx=newknt(sp);
   
    case 'optknt' %Gaffney/Powell, see matlab documentation
        %implement
    case 'chbpnt' %Chebychev-Demko points, see matlab documentation
        xx=chbpnt(sp.knots,sp.order);
    otherwise
        disp('Refinement method unknown');
        %error message/usage?
end

xx=setdiff(xx,origxx); %we keep only the new knots (syntax of fnrfn)
spr=fnrfn(sp,xx)  %xx is the new knots
xxstarnew=aveknt(spr.knots,spr.order);
N=length(xxstarnew);
mylegend=sprintf('Refined control Polygon, N=%d',N);
% waitforbuttonpress;
% pause(.3);
set(pr,'XData',xxstarnew,'YData',spr.coefs);%,'DisplayName',mylegend);
% legend({legsp,legcp,mystr,mylegend},'Position', [.7 .7 .2 .2]);
legend({legcp,mylegend},'Position', [.7 .7 .2 .2]);
% waitforbuttonpress;
% pause(.3);
set(pr,'Marker','o');
legend('off');
drawnow;

% xxuni=linspace(origxx(1),origxx(end),N);
% sampleValues=fnval(sp,xxuni); 
% this is simply sampling of the original spline as is with N points
% alternative: sample at the xxstar locations - compare performance?
splineValuesnew=fnval(spr,xxstarnew);
if nargin==7 && showSampledSpline
%     waitforbuttonpress;
%     pause(.3);
%     plot(xxuni,sampleValues,'m-o','DisplayName',mystr);
    mystr=sprintf('Sampling, N=%d',N);
    set(psample,'XData',xxstarnew,'YData',splineValuesnew,'DisplayName',mystr);
%     legend('show'); drawnow;
end
hold off;

end

