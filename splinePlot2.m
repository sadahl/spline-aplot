function spr =  splinePlot2(sp,insertFunction, plotSpline, noInsert)
%Plots the refined control polygon of sp using refinement method
%insertFunction

% Setting up the plot window
p=sp.order-1; % order of spline
scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/1.3 scrsz(3)/1.3 scrsz(4)/1.3])

%Plotting according to the dimensions of the spline coefficients
if sp.dim==1
    xxstar=aveknt(sp.knots,p+1); % initial sampling at the knot averages
    c1=xxstar;
    c2=sp.coefs;
    N=length(xxstar);
else % sp.dim==2
    c1=sp.coefs(1,:);
    c2=sp.coefs(2,:);
end % end if

legcp='Control Polygon';
pc = plot(c1, c2,'--ok');%, 'DisplayName', legcp); 
hold on;

spr=sp; % in case of no refinement

if nargin==2 || nargin==3 || (nargin>=4 && ~noInsert) %refinement
    %Choosing which knots to insert
    xx=insertFunction(sp); % example: xx = aveknt(sp.knots,sp.order)

    %Inserting knots
    spr=fnrfn(sp,xx) % spr holds the new spline after insertion of new knots
    xxstarnew=aveknt(spr.knots,spr.order);
    % N=length(xxstarnew);
    % legcpr=sprintf('Refined control Polygon, N=%d',N);
    % waitforbuttonpress;
    % pause(.3);
    if sp.dim==1
        c1r=xxstarnew;
        c2r=spr.coefs;
    else % sp.dim==2
        c1r=spr.coefs(1,:);
        c2r=spr.coefs(2,:);
    end %end if
    pr = plot(xxstarnew, spr.coefs,'-db');
    % set(pr,'XData',xxstarnew,'YData',spr.coefs);%,'DisplayName',mylegend);
    % legend({legsp,legcp,mystr,mylegend},'Position', [.7 .7 .2 .2]);
    % legend({legcp,mylegend},'Position', [.7 .7 .2 .2]);
    legend('off'); drawnow;

    if sp.dim==1
        knotPlotter(sp.knots, sp.coefs, xx); %not guaranteed to be 
            %the actually inserted knots due to fnrfn internal checks
    end
else %no refinement (no knots inserted)
   if sp.dim==1
        knotPlotter(sp.knots, sp.coefs);
   end
end % end if

% legsp='The spline';
if nargin>=3 && plotSpline
    fnpltpoints=fnplt(sp); %plotted with reference tool fnplt
    ps=plot(fnpltpoints(1,:),fnpltpoints(2,:), 'r-'); hold off;
end


