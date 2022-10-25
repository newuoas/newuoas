load ITERDATA;

fid=fopen('PROB','r');
PROB=textscan(fid,'%s');
PROB=PROB{1};
fclose(fid);

fid=fopen('SOLVER','r');
SOLVER=textscan(fid,'%s');
SOLVER=SOLVER{1};
fclose(fid);
if (KS ~=size(SOLVER))
    fprintf('The number of solvers does not match the file \"SOLVER\".\n');
    return;
end

logplot=1;
T = zeros(KP*KD,KS);

for s = 1:KS
    for i = 1:KP*KD
        p = min(floor(i/KD)+1,ceil(i/KD));
        if (mod(i,KD) == 0)
            d = KD;
        else
            d = mod(i,KD);
        end
        T(i,s) = meaniter(s,p,d);
    end
%    T(:,s) = reshape(maxiter(s,:,:),KP*KD,1);
%    T(:,s) = reshape(miniter(s,:,:),KP*KD,1);
%    T(:,s) = reshape(stditer(s,:,:),KP*KD,1);
%    T(:,s) = reshape(rstditer(s,:,:),KP*KD,1);
%    T(:,s) = reshape(oneiter(s,:,:),KP*KD,1);
end

name='mean';
perf(T, logplot, name, NaN,SOLVER);
data(T./(DIMENSION+1), name, NaN,SOLVER);

for s = 1:KS
    for i = 1:KP*KD
        p = min(floor(i/KD)+1,ceil(i/KD));
        if (mod(i,KD) == 0)
            d = KD;
        else
            d = mod(i,KD);
        end
         T(i,s) = stditer(s,p,d);
    end
end
name='std';
sperf(T, logplot, name, NaN,SOLVER);

for s = 1:KS
    for i = 1:KP*KD
        p = min(floor(i/KD)+1,ceil(i/KD));
        if (mod(i,KD) == 0)
            d = KD;
        else
            d = mod(i,KD);
        end
         T(i,s) = rstditer(s,p,d);
    end
end
name='rstd';
rsperf(T, logplot, name, NaN,SOLVER);



fid=fopen('mean','w');
for p = 1:KP
    for s = 1:KS
        for d = 1:KD
            fprintf(fid,'&%d',round(meaniter(s,p,d)));
        end
        fprintf(fid,'\\\\');
        fprintf(fid,'\n');
    end
end
fclose(fid);

fid=fopen('std','w');
for p = 1:KP
    for s = 1:KS
        for d = 1:KD
            fprintf(fid,'&%d',round(stditer(s,p,d)));
        end
        fprintf(fid,'\\\\');
        fprintf(fid,'\n');
    end
end
fclose(fid);

fid=fopen('rstd','w');
for p = 1:KP
    for s = 1:KS
        for d = 1:KD
            fprintf(fid,'&%3.2f',rstditer(s,p,d));
        end
        fprintf(fid,'\\\\');
        fprintf(fid,'\n');
    end
end
fclose(fid);

fid=fopen('mean-std','w');
for p = 1:KP
    for s = 1:KS
        for d = 1:KD
            fprintf(fid,'&%d\\,/\\,%d',round(meaniter(s,p,d)),round(stditer(s,p,d)));
        end
        fprintf(fid,'\\\\');
        fprintf(fid,'\n');
    end
end
fclose(fid);

fid=fopen('mean-rstd','w');
    fprintf(fid,'\\FloatBarrier\n\\newpage\n'); 
    for p = 1:KP
        if(mod(p,10)==1)
        fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid,'\n');
        if (p == 1)
        fprintf(fid,'\\begin{table}\n\\centering\n\\topcaption{\\newuoa 与 \\newuoas 的数值表现 (迭代自然终止)}\n\\label{tab:newuoas}\n\\begin{tabular}{lllll}\n\\toprule &50 & 100 & 150 & 200\\\\\n\\midrule\n');
        else
        fprintf(fid,'\\begin{table}\n\\centering\n{续表 \\ref{tab:newuoas}: \\newuoa 与 \\newuoas 的数值表现 (迭代自然终止)}\\\\[\\Bcaptionskip]\n\\begin{tabular}{lllll}\n\\toprule &50 & 100 & 150 & 200\\\\\n\\midrule\n');
        end
        end
        for s = 1:KS
            if (mod(s,2) == 1) 
                fprintf(fid,'\\multirow{2}{\\namewidth}{%s}',PROB{p});
            else
                fprintf(fid,'                                    ',PROB{p});
            end 
            for d = 1:KD
                fprintf(fid,'&%d\\,/\\,%3.2f',round(meaniter(s,p,d)),rstditer(s,p,d));
            end
            fprintf(fid,'\\\\');
            fprintf(fid,'\n');
            if (mod(s,2) == 0) 
                if (mod(p,10) == 0) 
                    fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table}\n');
                else
                    fprintf(fid,'[\\ski]\n');
                end
            end
        end
    end
    fclose(fid);

fid=fopen('std-rstd','w');
for p = 1:KP
    for s = 1:KS
        for d = 1:KD
            fprintf(fid,'&%d\\,/\\,%3.2f',round(stditer(s,p,d)),rstditer(s,p,d));
        end
        fprintf(fid,'\\\\');
        fprintf(fid,'\n');
    end
end
fclose(fid);

fid=fopen('fail','w');
for p = 1:KP
    for s = 1:KS
        for d = 1:KD
            fprintf(fid,'&%d',fail(s,p,d));
        end
        fprintf(fid,'\\\\');
        fprintf(fid,'\n');
    end
end
fclose(fid);

for s = 1:KS
    fprintf('Number of failures of %s:',SOLVER{s});
    reshape(fail(s,:,:),KP,KD)
end


clear;




