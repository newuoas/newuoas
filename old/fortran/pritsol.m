load PRITERDATA;
logplot = 1;
T = zeros(KP*KD,KS);
sol=[1,2];

fid=fopen('SOLVER','r');
SOLVER=textscan(fid,'%s');
SOLVER=SOLVER{1};
fclose(fid);
if (KS ~=size(SOLVER))
    fprintf('The number of solvers does not match the file \"SOLVER\".\n');
    return;
end

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
    perf(T(:,sol), logplot,name,tau(pr),SOLVER(sol));
    data(T(:,sol)./(DIMENSION(:,sol)+1), name,tau(pr),SOLVER(sol));
    
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
    perf(T(:,sol), logplot,name,tau(pr),SOLVER(sol));

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
    perf(T(:,sol), logplot,name,tau(pr),SOLVER(sol));
    
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
    for p = 1:KP
        for s = 1:KS
            for d = 1:KD
                fprintf(fid,'&%d\\,/\\,%3.2f',round(meaniter(pr,s,p,d)),rstditer(pr,s,p,d));
            end
            fprintf(fid,'\\\\');
            fprintf(fid,'\n');
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





