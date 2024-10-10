load ITERDATA;
logplot=1;
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
perf(T(:,sol), logplot, name, NaN,SOLVER(sol));
data(T(:,sol)./(DIMENSION(:,sol)+1), name, NaN,SOLVER(sol));

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
perf(T(:,sol), logplot, name, NaN,SOLVER(sol));

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
perf(T(:,sol), logplot, name, NaN,SOLVER(sol));



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
for p = 1:KP
    for s = 1:KS
        for d = 1:KD
            fprintf(fid,'&%d\\,/\\,%3.2f',round(meaniter(s,p,d)),rstditer(s,p,d));
        end
        fprintf(fid,'\\\\');
        fprintf(fid,'\n');
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




