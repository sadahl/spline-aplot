%% Using the MATLAB Curve Fitting Toolbox - examples
% Studying the effect of one knot insertion (nonuniform)
% Simple example: quadratic B-spline, two non-zero coefficients
% Study the control polygon refinement in this simple case

% setting up the spline
p=2; %degree
knots=[0 .25 1]; % knots without multiplicities
coeffs=[0 2 0 1]; % spline coefficients
sp=spmak(augknt(knots,p+1),coeffs)% spmak uses order=degree+1
                                  % augknt ensures (p+1)-regular k.vector

% setting up window                                  
scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/1.3 scrsz(3)/1.3 scrsz(4)/1.3])

% initial sampling at the knot averages
xxstar=aveknt(sp.knots,p+1) 
N=length(xxstar);
legcp='Control Polygon';
pc = plot(xxstar, sp.coefs,'--ok');%, 'DisplayName', legcp); 
hold on;

%Choosing which knots to insert
xx=[.125 .75];

%Inserting knots using built-in function fnrfn
spr=fnrfn(sp,xx) % spr holds the new spline after insertion of new knots
xxstarnew=aveknt(spr.knots,spr.order)
N=length(xxstarnew);
% legcpr=sprintf('Refined control Polygon, N=%d',N);
% waitforbuttonpress;
% pause(.3);
pr = plot(xxstarnew, spr.coefs,'-db');
% set(pr,'XData',xxstarnew,'YData',spr.coefs);%,'DisplayName',mylegend);
% legend({legsp,legcp,mystr,mylegend},'Position', [.7 .7 .2 .2]);
% legend({legcp,mylegend},'Position', [.7 .7 .2 .2]);
legend('off'); drawnow;

knotPlotter(sp.knots, sp.coefs, xx);

fnpltpoints=fnplt(sp); %plotted with reference tool fnplt
% legsp='The spline';
ps=plot(fnpltpoints(1,:),fnpltpoints(2,:), 'r.'); hold off;
% saveas(gcf,'../figures/quadratic_insert','epsc');
