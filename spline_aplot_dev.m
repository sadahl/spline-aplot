function [ spr, tadded, flag ] = spline_aplot_dev( sp, C, R, S, Sigma )
%[ spr, tadded, flag ] = spline_aplot( sp, C, R, N, S, Sigma )
%   Returns the refined spline spr whose control points represent an
%   adequate plot of the original spline sp. Both sp and spr are in B-form.
%   Returns also the vector of added knots, tadded. (unsorted)
%   flag indicates why the function returned:
%       -1 (error)
%       0 (max knots: number of knots inserted reached maximum allowed
%       1 (tau reached):  exit due to C criterion reached
%       2 (epsilon reached): exit due to refinement acc. to cpdiff reached
%
%   Criterion C (scalar): finds the region of the cpolygon assumed to hold
%   greatest potential for visual improvement. Searches through the control
%   points and returns an index mu and a range parameter r. Choices:
%       1: length_simple
%       2: lengthsq
%       3: length_ratio (incomplete - not fully tested)
%       4: angledot_simple (default)
%       5: baselineDist
%       6: angledot
%       7: baselineDist_extended    % to be implemented
%
%   Rule R (scalar): determines how to find the value of the knot z
%   to be inserted in the current iteration. Choices:
%       0:   Inserts knot at the previous knot average.
%       0.5: Inserts knot to produce new cp midway btwn c_mu and c_{mu+1}.
%            (Default)
%       1:   Inserts knot at the knot average.
%       w:   Inserts knot with this barycentric value btwn knot averages. 
%       To analyze: better to insert midway between knots or at knots?
%       Ref. even degree, odd degree.
%
%   Insertion I: inserts the knot. Relies on MATLAB built-in function fnrfn
%   implementing Bohms algorithm.
%
%   Stopping criterion S (array): threshold values linked to visual quality,
%   governs the iterations. S = [ N, delta, skip, consec, tbreak, tau ]
%       N      maximum allowed number of knots to be inserted.
%              0 means no limit. Default: 500.
%       delta  threshold for minimal change in control polygon before
%              break, measured in pixels. 
%              Linked to epsilon=delta/ro*deltac where ro is screen size in
%              pixels and deltac is size of box holding the spline sp.
%              Default: 3 pixels.
%       skip   tests for eps criteria every skip iterations. Default: 1.
%       consec number of consecutive iterations yielding change below eps
%              needed before break. Default: 5.
%       tbreak bool (def: true), activate break threshold for the relevant
%              property chosen in criterion C. Iterations will stop if no
%              region under the threshold is found. Uses default tau:
%                   length methods:     epsilon
%                   ratio method:       1 + .1*delta/l_max
%                       (takes l_max = ro if no prephase)
%                   angle methods:      2*pi/180 (rad, corr. 2 degrees)
%                   distance methods:   epsilon (comp. can be refined)
%       tau    (optional): overwrites the default values specified under
%              tbreak.
%
%   Partitioning Sigma of the control polygon: 
%       array   If provided, determines the control point indices
%               that serve as break points between the local control,               
%               including the end point indices (1 and n).
%       scalar  Regular partitioning with local control polygons of 
%               the specified size, starting from the left. Value 0
%               is interpreted as p (spline degree). Default: 3.
%       

% To implement: checks that input is as expected. (p+1)-regular knots and
% no knot with multiplicity higher than p-1.
% check_regular
% check_continuous

% To implement: reuse of candidates in myC, so we don't recompute all
% candidates in every iteration, but use knotmu (location of inserted knot)
% to update only the values that were changed.
% Insert then only the new knots into myC and keep track of correct
% indices.
% This is an optimization comment, it does not change the method/poc. 

% Initializing counters
Ni=0; % number of inserted knots
myconsec=0; % number of consecutive checks on epsilon succeeded
spr=sp;  %ensures proper exit
tadded=NaN; % ensures proper exit
flag = -1;


% Checking for the knot placement rule R:
if nargin <3, R=0.5; end

% Setting up stopping criterias, possibly provided by the user
if nargin>=4
    N = S(1);
    delta = S(2);
    skip = S(3);
    consec = S(4);
    tbreak = S(5);
    if length(S)==6
        tau = S(6);
        usertau = 1;
    else
        usertau = 0;
    end
else
    disp('Using default values for stopping criterion S...');
    N = 10;%500;
    delta = 3;
    skip = 1;
    consec = 5;
    tbreak = 1;
    % tau is specified according to the method/criterion C
    usertau = 0;
end

% Setting up visual quality variables
scrsz = get(groot,'ScreenSize'); % screensize: [left bottom width height].
ro=norm(scrsz(3:4));
%ro=max(scrsz(3:4))+1/2*min(scrsz(3:4)) %alternative (approx. of diagonal)
if sp.dim == 1
    cpoints=[aveknt(sp.knots,sp.order);sp.coefs]; %spline function
else
    cpoints=sp.coefs; %spline curve
end
ll=min([cpoints']);
ur=max([cpoints']);
deltac=norm(ur-ll); % size of box containing the spline to plot (euclidean)
epsilon=delta/ro*deltac; % tolerance converted from pixels to eucl. dist.
p=sp.order-1;

% Preparing the geometric criterion C
switch C
    case 1 % length_simple
        if ~usertau, tau = epsilon; end
        r=0; % range parameter 
        myC = @C_length_simple;
    case 2 % lengthsq
        if ~usertau, tau = epsilon^2; end
        r=0; % range parameter
        myC = @C_lengthsq;
    case 3 % length_ratio
        if ~usertau, tau = (1 + .1*delta/ro); end
        r=0;
        myC = @C_lengthratio;
    case 4 % angledot_simple
        if ~usertau
            % see thesis ChapterX for this value of tau
            tauconvert = 1.543;
            alphatau = 2; %degrees
            tau = tauconvert*(alphatau*pi/180)^2;
        end            
        r=0; % range parameter 
        myC = @C_angledot_simple;
    case 5 % baselineDist
        if ~usertau, tau = epsilon^2; end  %baselineDist returns dist sq.
        r=0;
        myC = @C_baselineDist;
    case 6 % angledot
        if ~usertau
            alphatau = 2; % degrees
            tau = alphatau*pi/180
        end
        r=0; % range parameter 
        myC = @C_angledot;
%     case 7 % baselineDist_extended % to be implemented
%         if ~usertau, tau = epsilon^2; end
%         r=0;
        %         myC = @C_;
    otherwise
        disp('Unknown method or method not yet implemented');
        disp('Using the simplified angle method, C=4');
        if ~usertau, tau = 2*pi/180; end
        r=0; % range parameter 
        %         myC = @C_;
        
end

% Preparing the rule R
% p=sp.order-1; % degree of spline
% R = rem(p,2)/2; % 0 for even degree, 0.5 for odd degree?
% An idea, but theory closer to these ratios referring to knot vector, not
% knot averages.

%Prephase: Insert knots to ensure each control polygon segment is shorter
%than a specified value, for instance ro/2^6.
% l_max = ... % redefine tau if C=3 (lenghratio-method)
%To be implemented.

%Set up heap here? To be implemented.
%Heap selection managed by C method - 
%keep the indices of the 10% worst regions in an array.
%Must then alter the output of C to provide multiple candidates.

%Sigma: to be implemented.

%Main work
knotmu = 0; % no knot inserted yet.
lsize = 3; % size of local control polygons. One inner point.
while Ni < N
    %Searching using criteria C
    [mu, r, crit] = myC(cpoints, lsize, knotmu);
%     crit
    
    if tbreak && crit<tau % measured criteria within threshold
        disp('Criteria tau reached, exiting')
        flag=1;
        return;         % possible improvement: change to another crit. C.
    end       
    
    % Applying the rule R for finding where to place a knot
    % range parameter r currently not in use
    z=sp.knots(mu) + R*(sp.knots(mu+p)-sp.knots(mu));
        %See Section 3.3.1 in thesis for this equation
    
    % Knot insertion I, giving the new control points
    spr=fnrfn(sp,z);
    knotmu=findMu(sp.knots,z); 
        % (knotmu may also be extracted from the knot insertion method).
        % The value kan be reused in the next iteration of the while-loop.
        % Only the control points with index knotmu-p+1, ..., knotmu have
        % changed since the last iteration.
        % myC must then keep track of several (all) candidates
        % This requires som restructuring of the code. ( Use instead of r?)
    Ni=Ni+1;
    tadded(Ni)=z;
    if sp.dim == 1
        tstar=aveknt(spr.knots,p+1);
        cpointsnew=[tstar;spr.coefs]; %spline function
    else
        cpointsnew=spr.coefs; %spline curve
    end 
    
    % Stopping criteria S: Check the change inflicted
    if rem(Ni,skip)==0
        % location of z in knot vector
        % (index of knot immediately to the left) 
        myeps=cpdiff(cpoints,cpointsnew,knotmu,p);
        %myeps=cpdiff(cpoints,cpointsnew); % also possible
        % improvement potential: extract this from knot insertion process.
        if myeps<epsilon
            myconsec = myconsec + 1;
        else
            myconsec=0;
        end
        if myconsec >= consec % myeps<eps a number of consec times in a row
            disp('Reached the level of refinement specified (epsilon)')
            flag=2;
            return; % spline is judged to be refined enough
        end  
    end % end if
    sp=spr;
    cpoints=cpointsnew;
    
    % Consider changing the method C or the partitioning Sigma / lcps
    % after a number of iterations.
    
end % end while
disp('Maximum number of knots allowed were inserted')
flag=0;

end

