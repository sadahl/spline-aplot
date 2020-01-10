%Show the convergence of the control polygon

% Setting up a spline
breaks = 0:10; k=4; p=k-1; %cubic
tau = augknt(breaks,k);
taustar = aveknt(tau,k);
sp = spapi(tau,taustar,sin(taustar)); % interpolation
% taustar2 = aveknt(fnbrk(sp,'knots'),fnbrk(sp,'order'))
cvec = fnbrk(sp,'coefs');

% Inserting more knots, in two rounds
t = midpointRefine(tau);
tstar = aveknt(t,k);
bvec = Oslo2(p, tau, cvec, t); %bvec and t are associated
tmark = midpointRefine(t);
tmarkstar = aveknt(tmark,k);
bmarkvec = Oslo2(p, t, bvec, tmark); %bmarkvec and tmark are associated

%plotting
figure;
axes('LineWidth',1,'FontSize',14,'FontName','Arial');

%method 1: uniform sampling of spline function
N=100;
x=linspace(tau(1), tau(end), N);
fx=fnval(sp,x);
plot(x,fx,'k', 'DisplayName', 'The spline');
%other method: use matlab native plot function directly
% fnplt(sp);  %slower! - catches breaks
hold on;

plot(taustar, cvec, 'ko-', 'DisplayName', 'Original control polygon');
plot(tstar, bvec, 'b-', 'DisplayName', 'First refinement');
plot(tmarkstar, bmarkvec, 'm-', 'DisplayName', 'Second refinement');
knotPlotter(tau,cvec);
%insert more knotPlotter commands: show new knots in other colour
% axis([0 3 -1.5 1.5]);
legend('show','Location', 'northeast');
xlabel('t');ylabel('y');
title('Convergence of the control polygon');
hold off;