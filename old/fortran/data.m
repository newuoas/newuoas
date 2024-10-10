function data(T,name,tau, solver)
%DATA    Data profiles
%
% Zaikun ZHANG, Feb 2012

delsame = 0;
fontsize = 16;
linewidth = 2;
penalty = 2;
colors  = [ 'b' 'k'  'r'  'm' 'c' 'g' 'y'];
lines   = {'-' '-.' '--' '-.' '-'};
markers = ['>' ' ' 'o'  ' '  'v' '^' 'o'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (delsame == 1)
    %% Delete the problems for which all the solvers performs the same.
    [np,ns] = size(T);
    ind = true(np,1);
    for p = 1:np
        if (sum(isnan(T(p,:))) == ns || (sum(isnan(T(p,:))) == 0 && max(T(p,:)) <= min(T(p,:))+0.01))
            ind(p) = false;
        end
    end
    T = T(ind,:);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[np,ns] = size(T);
if(np == 0)  % Prevent T from being empty.
    T=ones(1,ns);
    np =1;
end

% Replace all NaN's with twice the max_ratio and sort.
max_data = max(max(T));
T(isnan(T)) = penalty*max_data;
T = sort(T);

% Plot stair graphs with markers.

clf;

h=figure(1);
for s = 1: ns
 [xs,ys] = stairs(T(:,s),(1:np)/np);
 option = [char(lines(s)) colors(s)];
%  option = ['-' colors(s) markers(s)];
 plot(xs,ys,option,'Linewidth',linewidth);
 hold on;
end

% Axis properties are set so that failures are not shown,
% but with the max_ratio data points shown. This highlights
% the "flatline" effect.

if (isnan(tau) ~= 1)
    str=strcat('\tau = 10^{',int2str(log10(tau)),'}');
%    title(str,'fontsize',fontsize);
end

xlabel('\alpha','fontsize',fontsize)
ylabel('\delta_s(\alpha)','fontsize',fontsize)

xlabh = get(gca,'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [0 .0175 0])

set(gca,'FontSize',fontsize);
axis([ 0 1.1*max(max_data,0.01) 0 1 ]);



% Legends and title should be added.
legend(solver, 'Location', 'SouthEast');

%epsname = strcat(name,'.eps');
%print(h,'-depsc', epsname);
%unix(['epstopdf ',epsname]);

set(h,'visible','off');
name = strcat(name,'-data');
saveas(h,name,'fig');
saveas(h,name,'epsc2');
epsname = strcat(name,'.eps');
unix(['epstopdf ',epsname]);



