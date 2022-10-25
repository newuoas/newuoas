function perf(T,logplot,name,tau,solver)
%PERF    Performace profiles
%
% PERF(T,logplot)-- produces a performace profile as described in
%   Benchmarking optimization software with performance profiles,
%   E.D. Dolan and J.J. More', 
%   Mathematical Programming, 91 (2002), 201--213.
% Each column of the matrix T defines the performance data for a solver.
% Failures on a given problem are represented by a NaN.
% The optional argument logplot is used to produce a 
% log (base 2) performance plot.
%
% This function is based on the perl script of Liz Dolan.
%
% Jorge J. More', June 2004
delsame = 0;

fontsize = 16;
linewidth = 2;
penalty = 2;

colors  = [ 'b' 'k'  'r'  'm' 'c' 'g' 'y'];
lines   = {'-' '-.' '--' '-.' '-'};
markers = ['>' ' ' 'o'  ' '  'v' '^' 'o'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(delsame==1)
% Delete the problems for which all the solvers performs the same.
    [np,ns] = size(T);
    ind = true(np,1);
    if (strncmpi(name,'rstd',4))
        for p = 1:np
            if (sum(isnan(T(p,:))) == ns || (sum(isnan(T(p,:))) == 0 && max(T(p,:)) <= min(T(p,:))))
                ind(p) = false;
            end
        end
    else
        for p = 1:np
            if (sum(isnan(T(p,:))) == ns || (sum(isnan(T(p,:))) == 0 && max(T(p,:)) <= min(T(p,:))+1))
                ind(p) = false;
            end
        end
    end
    T = T(ind,:);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[np,ns] = size(T);
if(np == 0)  % Prevent T from being empty.
    T=ones(1,ns);
    np =1;
end
T(~isnan(T)) = max(1.0D-6,T(~isnan(T)));
% Make sure that T does not have zero values.
% Think about the case when all the solvers perform so 
% stable that rsdt=0. If we do not do this, all the ratios in
% r will be NaN, which stands for failure.

% Minimal performance per solver
minperf = min(T,[],2);

% Compute ratios and divide by smallest element in each row.
r = zeros(np,ns);
for p = 1: np
  r(p,:) = T(p,:)/minperf(p);
end


if (logplot) 
   r = log2(r);
end

max_ratio = max(1.0D-6,max(max(r)));

% Replace all NaN's with twice the max_ratio and sort.
r(isnan(r)) = penalty*max_ratio;
r = sort(r);

% Plot stair graphs with markers.

clf;

h=figure(1);
for s = 1: ns
 [xs,ys] = stairs(r(:,s),(1:np)/np);
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
    %title(str,'fontsize',fontsize);
end

xlabel('log_2(\alpha)','fontsize',fontsize)
ylabel('\sigma_s(\alpha)','fontsize',fontsize)

xlabh = get(gca,'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [0 .0175 0])

set(gca,'FontSize',fontsize);
axis([ 0 1.1*max(max_ratio,0.01) 0 1 ]);



% Legends and title should be added.
legend(solver,'Location', 'SouthEast');

%epsname = strcat(name,'.eps');
%print(h,'-depsc', epsname);
%unix(['epstopdf ',epsname]);

set(h,'visible','off');
name = strcat(name,'-perf');
saveas(h,name,'fig');
saveas(h,name,'epsc2');
epsname = strcat(name,'.eps');
unix(['epstopdf ',epsname]);



