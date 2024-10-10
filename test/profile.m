function T = profile(frec, fmin, tau, n, testfeature)
%
% This version is not inteded to be released. It is only for test.
%
% All rights reserved.
%
% ZHANG Zaikun, 08/08/2016
% Department of Applied Mathematics, The Hong Kong Polytechnic University

[np, ns, nr, maxfun] = size(frec);
M = maxfun;

T = NaN(np, ns, nr);

f0 = -Inf(np, nr);
for ip = 1:np
    for ir = 1:nr
        f0(ip,ir) = frec(ip, 1, ir, 1);
    end
end

for ip = 1:np
    for is = 1:ns
        for ir = 1:nr
            if (min(frec(ip, is, ir, 1:M)) <= tau*f0(ip,ir) + (1-tau)*fmin(ip)) % Do not change the "if .. else ..." order, because frec(ip, is, ir, 1:M) may be a vector of NaNs.
                T(ip, is, ir) = find(frec(ip, is, ir, 1:M) <= tau*f0(ip,ir) + (1-tau)*fmin(ip), 1, 'first');
            else
                T(ip, is, ir) = NaN;
            end
        end
    end
end

T = mean(T, 3);

sol = textread('solvers', '%s');
for is = 1:ns
    sol{is} = regexprep(sol{is}, '_4test', '');
    sol{is} = regexprep(sol{is}, 'newuoa', 'NEWUOA');
    if strcmpi(sol{is},'fminunc')
        sol{is} = '\texttt{fminunc}';
    end
    if strcmpi(sol{is},'newuoasv03c')
        sol{is} = 'Algorithm 2';
    end
end


delsame = 0;

fontsize = 16;
linewidth = 1.2;
penalty = 2;

colors  = {'k' 'b'  '#0072BD' 'm' 'c' 'g' 'y' 'b' 'k' 'r'};
lines   = {'-' '-.' ':' '--' '-' '-' '-.' '-' '-.' '-'};
markers = ['>' ' ' 'o'  ' '  'v' '^' 'o'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(delsame==1)
% Delete the problems for which all the solvers performs the same.
    mask = true(np,1);
    for ip = 1:np
        if (sum(isnan(T(ip,:))) == ns || (sum(isnan(T(ip,:))) == 0 && max(T(ip,:)) <= min(T(ip,:))+1))
            mask(ip) = false;
        end
    end
    T = T(mask,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[np, ns] = size(T);
if(np == 0)  % Prevent T from being empty.
    T=ones(1,ns);
    np =1;
end

Tmin = min(T, [], 2);

r = zeros(np, ns);
for ip = 1: np
  r(ip, :) = T(ip, :)/Tmin(ip);
end
r = log2(r);
max_ratio = max(1.0D-6,max(max(r)));
r(isnan(r)) = penalty*max_ratio;
r = sort(r);

clf;
hfig=figure(1);
for is = 1:ns
    [xs,ys] = stairs(r(:,is), (1:np)/np);
    plot(xs, ys, 'LineStyle', lines{is},  'Color', colors{is}, 'Linewidth',linewidth);
    hold on;
end
axis([ 0 1.1*max(max_ratio,0.01) 0 1 ]);

%title(sprintf('Performance Profile with tolerance $\\tau=10^{%d}$', int32(log10(tau))), 'interpreter', 'latex');

% Legends and title should be added.
if (ns > 3)
    legend(sol,'Location', 'southeast','Orientation','vertical', 'interpreter', 'latex', 'fontsize', fontsize);
else
    %legend(sol,'Location', 'northoutside','Orientation','horizontal');
    legend(sol,'Location', 'southeast','Orientation','vertical', 'interpreter', 'latex', 'fontsize', fontsize);
end

xlabel('$\log_2(\alpha), \quad \alpha = \mathrm{NF}/\mathrm{NF}_{\min}$', 'fontsize', fontsize, 'interpreter', 'latex');
ylabel('$\pi_s(\alpha)$', 'fontsize', fontsize, 'interpreter', 'latex');
xlabh = get(gca,'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [0 0.016 0])
set(gca,'FontSize',fontsize);

set(hfig,'visible','off');
figname = strcat(testfeature, '_', int2str(n), '_', int2str(int32(-log10(tau))),'_perf');
epsname = strcat(figname,'.eps');
saveas(hfig, epsname, 'epsc2');
system('mkdir -p results');
system(['epstopdf ',epsname]);
system(['mv ', epsname, ' ./results']);
system(['mv ', strcat(figname, '.pdf'), ' ./results']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = NaN(ns, M);
for is = 1:ns
    for i = 1:M
        D(is, i) = length(find(T(:, is) <= i))/np;
    end
    for i = M+1:1.2*M
        D(is, i) = D(is, M);
    end
end

clf;
hfig=figure(2);
for is = 1:ns
    [xs, ys] = stairs([1:1.2*M]/(n+1), D(is, :));
    plot(xs, ys, 'LineStyle', lines{is},  'Color', colors{is}, 'Linewidth',linewidth);
    hold on;
end
axis([ 0 1.1*M/(n+1) 0 1 ]);

title(sprintf('Data Profile with tolerance $\\tau=10^{%d}$', int32(log10(tau))), 'interpreter', 'latex');

% Legends and title should be added.
if (ns > 3)
    legend(sol,'Location', 'southeast','Orientation','vertical', 'interpreter', 'latex', 'fontsize', fontsize);
else
    %legend(sol,'Location', 'northoutside','Orientation','horizontal');
    legend(sol,'Location', 'southeast','Orientation','vertical', 'interpreter', 'latex', 'fontsize', fontsize);
end

xlabel('$\beta = \mathrm{NF}/(n+1)$', 'fontsize',fontsize, 'interpreter', 'latex');
ylabel('$\delta_s(\beta)$', 'fontsize', fontsize, 'interpreter', 'latex');
xlabh = get(gca, 'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [0 0.016 0])
set(gca, 'FontSize', fontsize);

set(hfig,'visible','off');
figname = strcat(testfeature, '_', int2str(n), '_', int2str(int32(-log10(tau))),'_data');
epsname = strcat(figname,'.eps');
saveas(hfig, epsname, 'epsc2');
system('mkdir -p results');
system(['epstopdf ',epsname]);
system(['mv ', epsname, ' ./results']);
system(['mv ', strcat(figname, '.pdf'), ' ./results']);
