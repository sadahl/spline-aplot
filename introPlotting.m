% Experimenting with splines
% Studying the convergence of the control polygon of a spline.
% See: Matlab documentation -> Curve Fitting Toolbox -> Splines

%% Motivation: the spline lies within the control polygon
%From the control polygon, can you guess the shape of the spline?
%Plotting the spline only at the knot averages, then refining the sampling
%at the midpoints in several rounds.
%
% Spline originally in 'pp'-form
% Spline is later converted to B-form
% Assuming positive y-values (for static axis in plot)
xval=0:8; %xsize=length(x);
yvec=[3 5 1 4 5 8 3 4 4 2 0];
cs=spline(xval,yvec); % cubic spline interpolation (default)
% How to obtain the knots? Use fn2fm to express spline in B-form.
sp=fn2fm(cs,'B-');
str1='A random spline curve obtained from interpolation and its control polygon';
[spr,~,Nfinal]=splinePlot(sp,30,3,'uniform',str1, 1,0); hold on;
% p0 = plot(xval,yvec,'ro','DisplayName','Interpolation points'); % the original interpolation points
% legend('show'); drawnow; hold off;
hold on;
knotPlotter(spr.knots,spr.coefs);
% saveas(gcf,'../figures/random_01.pdf'); hold off; %check format and scaling!


%% Setting up another test-spline.
% Spline in 'B'-form
% Using spmak: takes knots and vectors as arguments.
p=3; % order of spline
xvec=0:8;
yvec=[3 5 1 4 5 8 3 4 4 2 0];
sp=spmak(augknt(xvec,p+1),yvec) % aptknt returns an acceptable knot-vector
str1='A random spline curve and its control polygon';
str1=' ';
[spr,~,~]=splinePlot(sp,25,3,'uniform',str1,1,1);
hold on;
knotPlotter(spr.knots,spr.coefs); hold off;
legend.Position=[0.2 0.6 0.1 0.2];
% axis on;
% saveas(gcf,'../figures/testprog_random','epsc');


%% Example from the MATLAB manual
% A sine curve
breaks = 0:10; k=4; %cubic
t = augknt(breaks,k);
x = aveknt(t,k);
sp = spapi(t,x,sin(x));
% tstar = aveknt(fnbrk(sp,'knots'),fnbrk(sp,'order'));
% a = fnbrk(sp,'coefs');
str1='Spline interpolation of sine curve and its control polygon';
[~,~,Nfinal]=splinePlot(sp,100,3,'uniform',str1, 1);


%% Simple example: cubic B-spline, one non-zero coefficient
% Study the control polygon refinement in this simple case
p=2; % order of spline
sp=spmak(augknt([0:1],p+1),[0 1 0])
str1='Spline with only one active coefficient and its control polygon';
[spr,~,Nfinal]=splinePlot(sp,16,4,'uniform',str1,1,1);
hold on;
knotPlotter(spr.knots, spr.coefs); hold off;
% str2='Spline with only one active coefficient and its control polygon';
% [~,~,Nfinal]=splinePlot(sp,100,4,'newknt',str2, 1); %only stnd values - improve
% str3='Spline with only one active coefficient and its control polygon';
% [~,~,Nfinal]=splinePlot(sp,100,4,'chbpnt',str3, 1);


%% Combining spline functions using fncmb
% Linear combinations, scaling, 

%% Refining the spline using fnrfn (standard knot insertion algorithm)
k=4;
sp1=spmak(augknt([0:1],k),[0 1 1 0])
str1='A cubic spline...';
splinePlot(sp1,20,4,'uniform',str1,1,1);
% ccc = fnbrk(fnrfn(circ,.5:4),'c');
% plot( aveknt( fnbrk(sp3,'knots'),k), fnbrk(sp3,'coefs'), 'r') 

%% Studying convexity properties of a cubic spline
p=2; % order of spline
sp=spmak(augknt([0:3],p+1),[0 2 1 3 5]) % alt: aptknt returns an acceptable knot-vector
% sp=spmak(augknt([0:3],p+1),[0 0 1 0 0])
str1='Quadratic spline with a concave and a convex part';
[spr,~,Nfinal]=splinePlot(sp,20,4,'uniform',str1,1,1);
hold on;
knotPlotter(spr.knots, spr.coefs); hold off;

%% Multivariate splines - generating a random spline
% First variable: 11-7=4 (cubic), second variable: 9-6=3 (quadr.)
fnplt( spmak({augknt(0:4,4),augknt(0:4,3)}, rand(7,6)) );

%% Studying Chaikin's algorithm
% Testing chaikin.m
% Chaikin fails
% Convergence limit of chaikin: not the spline curve

%Setting up a spline
%original test-curve, dimension=1
p=2; % order of spline
sp=spmak(augknt([0 .25 1],p+1),[0 2 0 1]);
c1=aveknt(sp.knots,p+1)
c2=sp.coefs

%alternative test-curve, dimension=2
p=3;
sp=sptest1
c1=sptest1.coefs(1,:);
c2=sptest1.coefs(2,:);
% sp=spmak(augknt([0 1],p+1),[0 8 0 1])
% sp=spmak(augknt([0:1],p+1),[0 1 1 0])
% sp=spmak(augknt([0:3],p+1),[0 2 1 3 5])

%Setting up a plot window
scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/1.3 scrsz(3)/1.3 scrsz(4)/1.3])

%Plotting original control polygon
pc = plot(c1,c2,'--ok'); hold on;
% knotPlotter(sp.knots, sp.coefs);

%Plotting the spline 'exactly'
fnpltpoints=fnplt(sp); %plotted with reference tool fnplt
ps=plot(fnpltpoints(1,:),fnpltpoints(2,:), 'r-');

%Preparing Chaikin's method
N=10; %maximum number of rounds
delta = 10; % ten pixels as tolerance
cdelta=cpdiag(c1,c2)
eps=delta/max(scrsz)*cdelta
cdiff=1; %result of the test on new and old coordinates
counter=0; %counter for nb of rounds 

while counter < N && cdiff<eps
[c1n,c2n]=chaikin(c1,c2);
counter=counter+1;
cdiff=cpdiff(c1,c2,c1n,c2n);
c1=c1n; c2=c2n;
end

% hold on;
pr = plot(c1n,c2n,'-b'); hold off;
% saveas(gcf,'../figures/chaikinfail','epsc');
% saveas(gcf,'../figures/chaikinfail2','epsc');



%% Studying Chaikin's algorithm - counting operations
%Number of arithmetic operations
%Number of computed points

n=5; % nb of orig points
kappa=6; %nb of rounds
Q=zeros(1,kappa);
N=zeros(1,kappa);
for k=1:kappa
    Q(k)=12*k+12*(n-2)*(2^k-1)
    N(k)=(n-2)*2^k+2;
end %end for

hold off;
plot(1:kappa,Q./N);

%% Studying Lane & Riesenfeld algorithm
% Testing lanerisenfeld.m
% Lane-Riesenfeld fails
% Convergence limit of the algorithm: not the spline curve

%Setting up a spline
p=2; % order of spline
% sp=spmak(augknt([0 .25 1],p+1),[0 2 0 1]) %used for chaikin, p=2
sp=spmak(augknt([0 1],p+1),[0 3 1])
% sp=spmak(augknt([0:1],p+1),[0 1 1 0])
% sp=spmak(augknt([0:3],p+1),[0 2 1 3 5])
c1=aveknt(sp.knots,p+1)
c2=sp.coefs

%Setting up a plot window
scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/1.3 scrsz(3)/1.3 scrsz(4)/1.3])

%Plotting original control polygon
pc = plot(c1,c2,'--ok'); hold on;
knotPlotter(sp.knots, sp.coefs);

%Plotting the spline 'exactly'
fnpltpoints=fnplt(sp); %plotted with reference tool fnplt
ps=plot(fnpltpoints(1,:),fnpltpoints(2,:), 'r-');

%Preparing Chaikin's method
N=3; %maximum number of rounds
delta = 10; % ten pixels as tolerance
cdelta=cpdiag(c1,c2)
eps=delta/max(scrsz)*cdelta
% cdiff=1; %result of the test on new and old coordinates
counter=0; %counter for nb of rounds 

[c1,c2]=laneriesenfeld(c1,c2,1);

while counter < N-1 %&& cdiff<eps
[c1n,c2n]=laneriesenfeld(c1,c2);
counter=counter+1;
% cdiff=
c1=c1n; c2=c2n;
end

% hold on;
pr = plot(c1n,c2n,'-b'); hold off;
% saveas(gcf,'../figures/lanefail','epsc');

%% Studying the effect of one knot insertion (odd)
% Simple example: cubic B-spline, two non-zero coefficients
% Study the control polygon refinement in this simple case
p=3; % order of spline
sp=spmak(augknt([0:1],p+1),[0 2 1 2]) % aptknt returns an acceptable knot-vector
% sp.knots
% aveknt(sp.knots,p+1)
% str1='Spline with only one active coefficient and its control polygon';
scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/1.3 scrsz(3)/1.3 scrsz(4)/1.3])

xxstar=aveknt(sp.knots,p+1); % initial sampling at the knot averages
N=length(xxstar);
legcp='Control Polygon';
pc = plot(xxstar, sp.coefs,'--ok');%, 'DisplayName', legcp); 

origxx=sp.knots; % keeping the original knots
xx=origxx;
xx=.5;  %choosing where to insert a knot
% xx=xxstar(3);

% legcpr=sprintf('Refined control Polygon, N=%d',N);
hold on;
pr = plot(xxstar, sp.coefs,'-db');%, 'DisplayName', mylegend); 

% xx=setdiff(xx,origxx) %we keep only the new knots (syntax of fnrfn)
spr=fnrfn(sp,xx)  %xx is the new knots % spr holds the new spline after insertion of new knots
xxstarnew=aveknt(spr.knots,spr.order);
N=length(xxstarnew);
legcpr=sprintf('Refined control Polygon, N=%d',N);
% waitforbuttonpress;
% pause(.3);
set(pr,'XData',xxstarnew,'YData',spr.coefs);%,'DisplayName',mylegend);
% legend({legsp,legcp,mystr,mylegend},'Position', [.7 .7 .2 .2]);
% legend({legcp,mylegend},'Position', [.7 .7 .2 .2]);

% set(pr,'Marker','o');
legend('off');
drawnow;

knotPlotter(sp.knots, sp.coefs, xx);

fnpltpoints=fnplt(sp); %plotted with reference tool fnplt
legsp='The spline';
ps=plot(fnpltpoints(1,:),fnpltpoints(2,:), 'r'); hold off;
% saveas(gcf,'../figures/cpbehaveodd','epsc');


%% Studying the effect of one knot insertion (even)
% Simple example: cubic B-spline, two non-zero coefficients
% Study the control polygon refinement in this simple case
p=2; % order of spline
sp=spmak(augknt([0 .5 1],p+1),[0 2 0 1])
scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/1.3 scrsz(3)/1.3 scrsz(4)/1.3])

xxstar=aveknt(sp.knots,p+1) % initial sampling at the knot averages
N=length(xxstar);
legcp='Control Polygon';
pc = plot(xxstar, sp.coefs,'--ok');%, 'DisplayName', legcp); 
hold on;

%Choosing which knots to insert
xx=[.25 .75];

%Inserting knots
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
legsp='The spline';
ps=plot(fnpltpoints(1,:),fnpltpoints(2,:), 'r'); hold off;

% saveas(gcf,'../figures/cpbehaveeven','epsc');

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


%% Testing new splinePlot2 function
p=2; % order of spline
sp=spmak(augknt([0 .25 .3 1],p+1),[0 3 .5 0 1])
spr1=splinePlot2(sp,@knotAverages,1);
spr2=splinePlot2(sp,@knotInteriors,1);

%% Testing new splinePlot2 function
% Test curve 0 ?
p=3; % order of spline
sptest0=spmak(augknt([0 .25 .3 1],p+1),[0 3 .5 1 0 1])
% spr1=splinePlot2(sp,@knotAverages,1);
% spr2=splinePlot2(sp,@knotInteriors,1);
splinePlot2(sptest0,1,1,1);
% tspother=[sp.knots(1:6) .3 .3 .3]  %testing effect on zoom
% cspother=[sp.coefs(1:3)]
% spother=spmak(tspother,cspother)
% hold on; fnplt(spother);

%% Test curve 1
% From Eriksen, p53.
p=3; n=16;
t=augknt(1:14,p+1);
c=[ [1 3 4 5 8 11 13 9 12 12 11 12.5 5.5 0 4 10],
    [2 1 2 3 6  9 11 11 7  4  0    2 4.5 8 15 4]];
sptest1=spmak(t,c);
splinePlot2(sptest1,1,1,1); % hakkete!
% hold on; fnplt(sp);
axis off;
% saveas(gcf,'../figures/testcurve1','epsc'); hold off;

%% Function, high degree, many variations
% Lyche, From Spline Methods Draft, 2018
p=25; n=26;
t=augknt(linspace(0,1,2),p+1);
c=[2 1 5 5 2 5 6 6 0 5 2 3 2 2 2 5 6 3 1 0 5 5 1 6 2 2];
sptest2=spmak(t,c)
splinePlot2(sptest2,1,1,1);
axis off;
% saveas(gcf,'../figures/testcurve2','epsc'); hold off;

%% Experimenting with 'imitating Chaikin's algorithm' (Refinement rule)
% Using test splines 0 and 1
% We insert knots such that we obtain new control points
% at locations corresponding to halfway or approx. halfway
% on the control polygon edges.

%% Visualizing results on the condition number
k=3; %cubic
uk(1,.5,k);
s1=linspace(0,.5,10);
s2=linspace(.5,1,10);
figure;
if mod(k,2)==0
plot(s1,uk(k/2-1,s1,k));hold on;
plot(s2,uk(k/2-1,1-s2,k));
else
plot(s1,uk((k-1)/2,1-s1,k));hold on;
plot(s2,uk((k-1)/2,s2,k));
end

%% Uniform sampling - sampling, waste of resources
p=4; % order of spline
sp=spmak(augknt([0  .1 .2 .3 1],p+1),[0 1 1 0 0 0 0])
% spr1=splinePlot2(sp,@knotAverages,1);
% spr2=splinePlot2(sp,@knotInteriors,1);
splinePlot2(sp,1,1,1);
% saveas(gcf,'../figures/samplewaste','epsc'); hold off;

%% Long arc length - uniform sampling inadequate. 'Flower'.
p=2; n=25;
t=augknt(1:n-p+1,p+1);
% figure();
% axis([0 15 0 15]);
% [x,y]=ginput(n); %get new points
% c=[x,y];
% c=c';

c=[[7.5 5 0 7.5 0 0 7.5 0 5 7.5 5 10 7.5 10 15 7.5 15 15 7.5 15 10 7.5 10 5 7.5],
[7.5 0 5 7.5 5 10 7.5 10 15 7.5 15 15 7.5 15 10 7.5 10 5 7.5 5 0 7.5 0 0 7.5]];

sptestlong=spmak(t,c)
splinePlot2(sptestlong,1,1,1);
% splinePlot2(sptest1,1,1,1);
hold on;
N=1000;
% tt=linspace(min(sptest1.knots),max(sptest1.knots),N);
tt=linspace(min(sptestlong.knots),max(sptestlong.knots),N);
finercurve=fnval(sptestlong,tt);
% finercurve=fnval(sptest1,tt);
sel=finercurve;
plot(sel(1,:),sel(2,:),'r-');
% saveas(gcf,'../figures/matlabntoosmall','epsc'); hold off;
% saveas(gcf,'../figures/matlabntoosmallflower','epsc'); hold off;

%% Complicated knot vector, simple spline / geometric shape
p=2; % order of spline
sp=spmak(augknt([0:1],p+1),[0 1 0])
insert = @(x) [0.1 .2];
spr=splinePlot2(sp,insert,1,0);
% saveas(gcf,'../figures/simplevscomplicated','epsc'); hold off;

%% Arccos and sqrt vs. other transformation
xx=linspace(-1,1,100); figure();hold on;
plot(xx,acos(xx),'r')
% plot(xx,-xx,'k--')
% plot(xx,-xx.^3./abs(xx),'k')
plot(xx,pi/2*-xx.*abs(xx)+pi./2,'b');
% saveas(gcf,'../figures/cosvariation','epsc'); hold off;
% myf=@(alpha) sin(pi/2-alpha); hold on;
% plot(xx,myf(xx),'m');

%% Sensibility of the angle transformation for the relevant limit values
% We must confirm the variations for alpha in [0, 10] (degrees), 
% since values in this interval are likely to be chosen as angle limits.
% It turns out a quadratic polynomial gives a good description of their
% relation.
alphatau=10; % degrees
alpha=linspace(0, alphatau*pi/180, 11);
vals=cos(alpha);
talpha=-pi/2*vals.*abs(vals)+pi/2;
figure();
plot(alpha*180/pi, talpha*180/pi, 'r'); % plot with axes in degrees
    % should fit well with a regular quadratic polynomial in this range
[f, gof] = fit(alpha', talpha', 'poly2')
% f = 
% 
%      Linear model Poly2:
%      f(x) = p1*x^2 + p2*x + p3
%      Coefficients (with 95% confidence bounds):
%        p1 =       1.543  (1.536, 1.55)
%        p2 =    0.002475  (0.001209, 0.003741)
%        p3 =  -3.129e-05  (-7.878e-05, 1.62e-05)
%   rsquare: 1.000
tauconvert = 1.543;
ftilde = @(x) tauconvert*x.^2;  % NB: x in radians
ftilde(alpha)
hold on;
plot(alpha*180/pi, ftilde(alpha)*180/pi, '--k'); % a good fit
% saveas(gcf,'../figures/angle_trnsflimits','epsc'); hold off;


%% Visual smoothness - curvature and arc length
% From Eriksen, p41
p=3;
t=[0 0 0 0 1 1 1 1];
c1=[[2.5 0.5 4.5 2.5],
    [1.5 2.1 2.1 1.5]];
add=zeros(2,4);
add(2,:)=1
% c1=c1+add;
c2=[[1 2 3 4],
    [1 1.1 1.1 1]];
% c2=c2-add;
sp1=spmak(t,c1);
sp2=spmak(t,c2);

%evaluating in only N points
N=15;
tt=linspace(0,1,N);
% points1=fnplt(sp1);
% points2=fnplt(sp2);
points1=fnval(sp1,tt);
points2=fnval(sp2,tt);

figure();
plot(points1(1,:),points1(2,:),'k'); hold on;
plot(points2(1,:),points2(2,:),'k');
plot([1 1.003],[2 2],'w'); %for scaling/output - prevents trimming
% axis([1 4 1 2]);
axis off;
% saveas(gcf,'../figures/vissmooth','epsc'); hold off;

% studying arc length for same curvature
c3=[[3 0 5 2],
    [1 2 2 1]];
c4=[[6 0 10 4],
    [1 3 3 1]];
c5=[[11 -1 19 7],
    [1 5 5 1]];
sp3=spmak(t,c3);
sp4=spmak(t,c4);
sp5=spmak(t,c5);

%evaluating in only NN points
NN=20;
tt=linspace(0,1,NN);
% points1=fnplt(sp1);
% points2=fnplt(sp2);
points3=fnval(sp3,tt);
points4=fnval(sp4,tt);
points5=fnval(sp5,tt);

figure();
plot(points3(1,:),points3(2,:),'k'); hold on;
plot(points4(1,:),points4(2,:),'k');
plot(points5(1,:),points5(2,:),'k');
% axis([2 11 1 4]);
axis off;
% saveas(gcf,'../figures/vissmoothscale','epsc'); hold off;

%% Studying approximations. Oscillations
p=3;
t=[0 0 0 0 .5 1 1 1 1];
x=[0 2 5 8 10];
y=[0 2 3 -2 1.6];
c1=[x;y]

sp1=spmak(t,c1)
points1=fnplt(sp1)
figure();
plot(points1(1,:),points1(2,:),'r');
axis off; hold on;

%extracting only NN points
N=length(points1);
NN=12;
skip=floor(N/NN);
sinvec=zeros(2,NN+1);
for i=0:NN-1
    sinvec(:,i+1)=points1(:,1+skip*i)
end
sinvec(:,NN+1)=points1(:,N);
osc=linspace(0,6*pi,NN+1)
sinvec(2,:)=sinvec(2,:)+.2*sin(osc)
plot(sinvec(1,:),sinvec(2,:),'k'); hold off;
% saveas(gcf,'../figures/oscillation','epsc'); hold off;

figure();
plot(points1(1,:),points1(2,:),'r'); hold on; axis off;
shiftvec=zeros(2,NN+1);
for i=0:NN-1
    shiftvec(:,i+1)=points1(:,1+skip*i)
end
shiftvec(:,NN+1)=points1(:,N);
% add=linspace(0,6*pi,NN+1);
halfway=floor(NN/2)
shiftvec(2,1:halfway)=shiftvec(2,1:halfway)+.15;
shiftvec(2,halfway+1:NN+1)=shiftvec(2,halfway+1:NN+1)-.15;
plot(shiftvec(1,:),shiftvec(2,:),'k');
% saveas(gcf,'../figures/offset','epsc'); hold off;

%% control polygon - mimics the shape
p=4;
coeffs=[3 5 1 4 5 8 3 4 4 2 0];
n=length(coeffs);
t=linspace(0,1,n-p+1)
sp=spmak(augknt(t,p+1),coeffs)
splinePlot2(sp,1,1,1);
% saveas(gcf,'../figures/cpsmooth','epsc'); hold off;

%% operation count for overall adaptive algorithm
N0vals=[10 100]; %nb of points to start with
% K=[10 50 100 200 500]; %nb of inserted knots
K=[10 50 400];% 100 500];
Nk=length(K);
p=[0:30]; %degree of spline
Q=53*p; % number of operations per iteration
figure();
N0=N0vals(1);
for i = 1:Nk
    R=Q*K(i)/(N0+K(i));
    plot(p,R,'-.k'); hold on;
end %end for
N0=N0vals(2);
for i = 1:Nk
    R=Q*K(i)/(N0+K(i));
    plot(p,R,'-.b'); hold on;
end %end for
Reval=3/2*p.*(p+1); % cost per point for the unif. eval. method, per dim.
% plot(p,Reval,'m'); % one dimension - not relevant here
plot(p,2*Reval,'r');
xlabel('Degree, p'); ylabel('Operations per final point, R')
axis([0 30 0 1500]);

% saveas(gcf,'../figures/opscompare','epsc'); hold off;

%% general operation count for overall adaptive algorithm
yvec=[.2 1 2 10 1000] % proportion of inserted knots rel. to nb of orig. kn.
p=[0:30]; %degree of spline
sveccolor=['-.k' '- k']  %keep three characters for each
figure();
for s = 2:2 %dimension of coordinates
    for i = 1:length(yvec)
        y=yvec(i);
    %     k=3/2*(p+1)*(y+1)*s/y-13-3*s
        k=3/2*(p+1)*s*(1+1/y)-3*s-13;
        plot(p,k,sveccolor((3*s-2):3*s)); hold on;
    end %end for
end % end for

%horizontal sections
l1=2
l2=5
a1=23
a2=40
b1=30
b2=56
lb=19
horiz=[l1 l2 lb a1 b1 a2 b2];
for i=1:length(horiz)
    plot([p(1) p(length(p))],[horiz(i) horiz(i)],'-.')
end % end for
    
xlabel('Degree, p'); ylabel('Threshold, k')
axis([0 30 0 70]);

% saveas(gcf,'../figures/opscomparegen','epsc'); hold off;

%% thresholds, y -> inf
l1=2
l2=5
lb=19
a1=23
a2=40
b1=30
b2=56
horiz=[l1 l2 lb a1 b1 a2 b2];
klim_s1=(horiz+29/2)*2/3
klim_s2=(horiz+16)/3
