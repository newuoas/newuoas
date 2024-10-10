KS = textread('SOLVERNUM'); % Number of solvers.
KD = textread('DIMNUM'); % Number of different dimensions.
KR = textread('RANDNUM'); % Number of random tests.

clear FEVL FOPT;
for i = 1:KS
    FEVL(i,:) = textread(strcat('N',int2str(i)));
    FOPT(i,:) = textread(strcat('F',int2str(i)));
end
FSTART=textread('FSTART');
DIM = textread('DIM');

fid=fopen('SOLVER','r');
SOLVER=textscan(fid,'%s');
SOLVER=SOLVER{1};
fclose(fid);
if (KS ~= length(SOLVER))
    fprintf('The number of solvers does not match the file \"SOLVER\".\n');
    return;
end

NTEST = size(FEVL,2); 
KP = NTEST/(KD*KR); % Number of problems.
k = 1;

%NPENALTY = 1; % Penalty on failure. 
%NPENALTY = 4; % Penalty on failure. 
NPENALTY = NaN; % Penalty on failure. 



iter = zeros(KS,KR);
fun = zeros(KS,KR);
fail = -Inf*ones(KS,KP,KD);
meaniter =  -Inf*ones(KS,KP,KD);
miniter =  -Inf*ones(KS,KP,KD);
maxiter =  -Inf*ones(KS,KP,KD);
stditer =  -Inf*ones(KS,KP,KD);
rstditer =  -Inf*ones(KS,KP,KD);
oneiter =  -Inf*ones(KS,KP,KD);
minfun = Inf*ones(KP,KD);
startfun = Inf*ones(KP,KD);


%tol = 1.0D-10;
%tol = 1.0D-8;
tol = 1.0D-6

for p = 1:KP
    for d = 1:KD
        for s = 1:KS
            iter(s,:) = FEVL(s,(0:KR-1)*KD+d+(p-1)*KD*KR);
            fun(s,:) = FOPT(s,(0:KR-1)*KD+d+(p-1)*KD*KR);
        end
        itermax = max(max(iter));
        minfun(p,d) = min(min(fun));
        fstart = FSTART((0:KR-1)*KD+d+(p-1)*KD*KR);
        startfun(p,d) = mean(fstart);
        fail(:,p,d) = 0;
        for r= 1:KR
            for s= 1:KS
                if ((fun(s,r)-minfun(p,d))/max(1,abs(minfun(p,d)))>=tol || 1-(startfun(p,d)-fun(s,r))/(startfun(p,d)-minfun(p,d))>=tol || minfun(p,d)>=startfun(p,d)*(1-tol))
                    iter(s,r) = itermax*NPENALTY;
                    fail(s,p,d) = fail(s,p,d)+1;
                end
            end
        end

        for s = 1:KS
            meaniter(s,p,d) = mean(iter(s,:));
            miniter(s,p,d) = min(iter(s,:));
            maxiter(s,p,d) = max(iter(s,:));
            stditer(s,p,d) = std(iter(s,:));
            rstditer(s,p,d) = stditer(s,p,d)./meaniter(s,p,d);
            oneiter(s,p,d) = iter(s,k);
        end
    end
end

DIMENSION=zeros(KP*KD,KS);
for i = 1:KS
    DIMENSION(:,i)=DIM;
end

save ITERDATA  KS KP KD KR DIMENSION  meaniter miniter maxiter stditer rstditer oneiter fail;




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
%    T(:,s) = reshape(maxiter(s,:,:),KP*KD,1);
%    T(:,s) = reshape(miniter(s,:,:),KP*KD,1);
%    T(:,s) = reshape(stditer(s,:,:),KP*KD,1);
%         T(i,s) = rstditer(s,p,d);
%    T(:,s) = reshape(oneiter(s,:,:),KP*KD,1);
    end
end

name='mean';
perf(T, logplot, name, NaN, SOLVER);
data(T./(DIMENSION+1), name, NaN, SOLVER);

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
sperf(T, logplot, name, NaN, SOLVER);

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
rsperf(T, logplot, name, NaN, SOLVER);

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





