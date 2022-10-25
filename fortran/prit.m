load PRITERDATA;

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

logplot = 1;
T = zeros(KP*KD,KS);

for pr = 1:KPR
    for s = 1:KS
        for i = 1:KP*KD
            p = min(floor(i/KD)+1,ceil(i/KD));
            if (mod(i,KD) == 0)
                d = KD;
            else
                d = mod(i,KD);
            end
            T(i,s) = meaniter(pr,s,p,d);
        end
    %    T(:,s) = reshape(maxiter(pr,s,:,:),KP*KD,1);
    %    T(:,s) = reshape(miniter(pr,s,:,:),KP*KD,1);
    %    T(:,s) = reshape(stditer(pr,s,:,:),KP*KD,1);
    %    T(:,s) = reshape(rstditer(pr,s,:,:),KP*KD,1);
    end
    name=int2str(pr);
    name=strcat('mean',name);
    perf(T, logplot,name,tau(pr),SOLVER);
    data(T./(DIMENSION+1), name,tau(pr),SOLVER);
    
    for s = 1:KS
        for i = 1:KP*KD
            p = min(floor(i/KD)+1,ceil(i/KD));
            if (mod(i,KD) == 0)
                d = KD;
            else
                d = mod(i,KD);
            end
            T(i,s) = rstditer(pr,s,p,d);
        end
    end

    name=int2str(pr);
    name=strcat('rstd',name);
    rsperf(T, logplot,name,tau(pr),SOLVER);

    for s = 1:KS
        for i = 1:KP*KD
            p = min(floor(i/KD)+1,ceil(i/KD));
            if (mod(i,KD) == 0)
                d = KD;
            else
                d = mod(i,KD);
            end
            T(i,s) = stditer(pr,s,p,d);
        end
    end

    name=int2str(pr);
    name=strcat('std',name);
    sperf(T, logplot,name,tau(pr),SOLVER);
    
    
    file=strcat('mean',int2str(pr));
    fid=fopen(file,'w');
    for p = 1:KP
        for s = 1:KS
            for d = 1:KD
                fprintf(fid,'&%d',round(meaniter(pr,s,p,d)));
            end
            fprintf(fid,'\\\\');
            fprintf(fid,'\n');
        end
    end
    fclose(fid);
    
    file=strcat('std',int2str(pr));
    fid=fopen(file,'w');
    for p = 1:KP
        for s = 1:KS
            for d = 1:KD
                fprintf(fid,'&%d',round(stditer(pr,s,p,d)));
            end
            fprintf(fid,'\\\\');
            fprintf(fid,'\n');
        end
    end
    fclose(fid);
    
    file=strcat('rstd',int2str(pr));
    fid=fopen(file,'w');
    for p = 1:KP
        for s = 1:KS
            for d = 1:KD
                fprintf(fid,'&%3.2f',rstditer(pr,s,p,d));
            end
            fprintf(fid,'\\\\');
            fprintf(fid,'\n');
        end
    end
    fclose(fid);
    
    file=strcat('mean-std',int2str(pr));
    fid=fopen(file,'w');
    for p = 1:KP
        for s = 1:KS
            for d = 1:KD
                fprintf(fid,'&%d\\,/\\,%d',round(meaniter(pr,s,p,d)),round(stditer(pr,s,p,d)));
            end
            fprintf(fid,'\\\\');
            fprintf(fid,'\n');
        end
    end
    fclose(fid);
    
    file=strcat('mean-rstd',int2str(pr));
    fid=fopen(file,'w');
    fprintf(fid,'\\FloatBarrier\n\\newpage\n'); 
    for p = 1:KP
        if(mod(p,10)==1)
        fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid,'\n');
        if (p == 1)
        fprintf(fid,'\\begin{table}\n\\centering\n\\topcaption{\\newuoa 与 \\newuoas 的数值表现 (精度 $\\tau = 10^{-%d}$)}\n\\label{tab:newuoas%d}\n\\begin{tabular}{lllll}\n\\toprule &50 & 100 & 150 & 200\\\\\n\\midrule\n',pr,pr);
        else
        fprintf(fid,'\\begin{table}\n\\centering\n{续表 \\ref{tab:newuoas%d}: \\newuoa 与 \\newuoas 的数值表现 (精度 $\\tau = 10^{-%d}$)}\\\\[\\Bcaptionskip]\n\\begin{tabular}{lllll}\n\\toprule &50 & 100 & 150 & 200\\\\\n\\midrule\n',pr,pr);
        end
        end
        for s = 1:KS
            if (mod(s,2) == 1) 
                fprintf(fid,'\\multirow{2}{\\namewidth}{%s}',PROB{p});
            else
                fprintf(fid,'                                    ',PROB{p});
            end 
            for d = 1:KD
                fprintf(fid,'&%d\\,/\\,%3.2f',round(meaniter(pr,s,p,d)),rstditer(pr,s,p,d));
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
    
    file=strcat('std-rstd',int2str(pr));
    fid=fopen(file,'w');
    for p = 1:KP
        for s = 1:KS
            for d = 1:KD
                fprintf(fid,'&%d\\,/\\,%3.2f',round(stditer(pr,s,p,d)),rstditer(pr,s,p,d));
            end
            fprintf(fid,'\\\\');
            fprintf(fid,'\n');
        end
    end
    fclose(fid);
    
    file=strcat('fail',int2str(pr));
    fid=fopen(file,'w');
    for p = 1:KP
        for s = 1:KS
            for d = 1:KD
                fprintf(fid,'&%d',fail(pr,s,p,d));
            end
            fprintf(fid,'\\\\');
            fprintf(fid,'\n');
        end
    end
    fclose(fid);
end


clear;





