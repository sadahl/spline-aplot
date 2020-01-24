% Defines test curves for experimenting with adaptive plotting.
% Specially constructed test splines are named: sptest0, ..., sptest4.
% Generates random splines, held in array randomsplines.

%% Simple quadratic Bezier function
p=2; n=3;
t=[0 0 0 1 1 1];
c=[0 1 1];
sptest0=spmak(t,c);
% splinePlot2(sptest0,1,1,1);

%% Function, quintic, many variations
% Lyche, From Spline Methods Draft, 2018
p=5; n=10;
t=augknt(linspace(0,1,6),p+1);
c=[5 0 10 6 11 8 6 3 12 8];
sptest1=spmak(t,c);
% splinePlot2(sptest1,1,1,1);

%% Function, high degree, many variations
% Lyche, From Spline Methods Draft, 2018
p=25; n=26;
t=augknt(linspace(0,1,2),p+1);
c=[2 1 5 5 2 5 6 6 0 5 2 3 2 2 2 5 6 3 1 0 5 5 1 6 2 2];
sptest2=spmak(t,c);
% splinePlot2(sptest2,1,1,1);

%% Curve, low degree, long arc, straight sections and curved sections
% From Eriksen, Plotting of splines, 1993
p=3; n=16;
t=augknt(1:14,p+1);
c=[ [1 3 4 5 8 11 13 9 12 12 11 12.5 5.5 0 4 10],
    [2 1 2 3 6  9 11 11 7  4  0    2 4.5 8 15 4]];
sptest3=spmak(t,c);
% splinePlot2(sptest3,1,1,1); % hakkete!
% hold on; fnplt(sp);
%axis off;

%% Curve, containing a small curl
% From Eriksen, Plotting of splines, 1993
% NB: contains succeding control points that are equal.
p=3; n=18;
t=[0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5];
c=[[4 5 1 2 2 3 3.998 4 4 4.002 3.999 4.001 4.001 4.003 6 7 7 8 6 5],
    [2 8 0.8 1 1 1.2 0.999 1 1 1.001 1.001 1 1 0.999 1.5 2 2 2.5 3 1]];
sptest4=spmak(t,c);
% splinePlot2(sptest4,1,1,1);

%% Curve, high degree?
% ... to be decided

%% Generation of random curves - 'havoc'
% Repeat NN times, for instance NN=1000:
NN=1000;
smax=2;
clear randomsplines; % First row: s=1. Second row: s=2...

for s=1:smax
    for i=1:NN
        % Generate random order n in range typically [3, 500] (integer)
        nlower=3; nupper=500;
        n = random('unid',nupper-nlower)+nlower-1;

        % Generate random degree p in range typically [2, 30] (integer)
        plower=2; pupper=30;
        p = min(random('unid', pupper-plower)+plower-1,n-1);

        % Generate random (p+1)-regular knot vector t in range 
        % typically [0 10], of length n+p+1
        tlower=0; tupper=10;
%         t = sort( random('unif', tlower, tupper, 1, n+p+1-2*p) )
        t = sort( random('unif', tlower, tupper, 1, n+p+1-2*p) );
%         t = [t(1)*ones(1,p) t t(n+p+1-2*p)*ones(1,p)];
        t = augknt(t, p+1); % augments the knot vector.
            % NB: does not check for inner multiplicities, which may be p.
            % (but not likely).
            % Possible improvement: generate multiplicities to ensure the
            % adaptive methods are tested against many such cases as well.

        % Generate random c in range typically [0, 10] / [0 0; 10 10]
        % We could also here opt for a normal distribution, but keep the
        % uniform distribution as we expect harder curves this way.
        clower=0; cupper=10;
        c = random('unif', clower, cupper, s, n);
        
        % Make the spline object
        randomsplines(s,i)=spmak(t, c);
    end
end

%% Generation of additional test curves - 'interpolated'
% Choose a pattern
% Interpolate using spapi function or similar
% ...