KS = textread('SOLVERNUM'); % Number of solvers.
KD = textread('DIMNUM'); % Number of different dimensions.
KR = textread('RANDNUM'); % Number of random tests.
DIM = textread('DIM');

clear FOPT;
for i = 1:KS
    FOPT(i,:) = textread(strcat('F',int2str(i)));
end
k = 1;

%NPENALTY = 4; % Penalty on failure. 
NPENALTY = NaN; % Penalty on failure. 

tau = 1.0D-1.^(1:1:10);
tau = [tau, NaN];
KPR = length(tau); % Number of precisions.
logplot=1;

fid=fopen('SOLVER','r');
SOLVER=textscan(fid,'%s');
SOLVER=SOLVER{1};
fclose(fid);
if (KS ~=size(SOLVER))
    fprintf('The number of solvers does not match the file \"SOLVER\".\n');
    return;
end

fid=fopen('PROB','r');
PROB=textscan(fid,'%s');
PROB=PROB{1};
fclose(fid);
KP =length(PROB); % Number of problems.
SUFFIX = '.FOPTREC';

colnum=zeros(KS,1);
FEVL=zeros(KS,KD*KR);
fail=zeros(KPR,KS,KP,KD);
iter=zeros(KPR,KS,KP,KD,KR);
meaniter=zeros(KPR,KS,KP,KD);
maxiter=zeros(KPR,KS,KP,KD);
miniter=zeros(KPR,KS,KP,KD);
stditer=zeros(KPR,KS,KP,KD);
rstditer=zeros(KPR,KS,KP,KD);
minf=zeros(KD*KR,1);
terminatef=zeros(KD*KR,1);


for p = 1:KP
    fprintf('\n Processing data for Problem %d: %s...', p, PROB{p});
    clear FRECTEMP;
    for s = 1:KS
        clear temp;
        temp=textread((strcat(PROB{p},SUFFIX,int2str(s))),'','emptyvalue',NaN);
        colnum(s) = size(temp,2);
        FRECTEMP(s,1:KD*KR,1:colnum(s))=temp;
    end
    FREC = NaN*ones(size(FRECTEMP));
    % size(FREC)=[KS,KD*KR,max(Function evaluations of the solvers for different
    % variants of Problem p)]. For each s (1<=s<=KS) and i (1<=i<=KP*KR), 
    % FREC(s,i,:) is the record of function values when solver s solves the i-th
    % variant of problem p, NaN following the function values so that FREC makes
    % an array.
    % 
    for s = 1:KS
        FREC(s,1:KD*KR,1:colnum(s))=FRECTEMP(s,1:KD*KR,1:colnum(s));
        for i = 1:KD*KR
            FEVL(s,i) = find(isnan(FREC(s,i,:)) ~= 1,1,'last');
        end
%        if (~isempty(find(FEVL(s,:)~=NF(s,(p-1)*KD*KR+(1:KD*KR)))))
%            fprintf('ERROR!');
%            FEVL(s,:) - NF(s,(p-1)*KD*KR+1:KD*KR)
%        end
    end
    for pr = 1:KPR-1
        for s = 1:KS
            for i = 1:KD*KR
                if (mod(i,KD) == 0)
                    d = KD;
                else
                    d = mod(i,KD);
                end
                r = min(floor(i/KD)+1,ceil(i/KD));
                minf(i)=min(min(min(FREC(:,d+(0:KR-1)*KD,:))));
                temp=find(FREC(s,i,:) <= tau(pr)*FREC(s,i,1)+(1-tau(pr))*minf(i),1,'first'); 
                if (isempty(temp) || minf(i)>=FREC(s,i,1)*(1-1.0D-6))
                    iter(pr,s,p,d,r)=NPENALTY*max(max(FEVL(:,d+(0:KR-1)*KD)));
                    fail(pr,s,p,d) = fail(pr,s,p,d)+1;
%                    fprintf('FAILURE: %f,%d,%d,%d\n',tau(pr),s,p,i);
                else
                    iter(pr,s,p,d,r)=temp;
                end
            end
        end
    end
    for pr = KPR  % For each problem, let the worst optimal value found be the
                  % threshold for termination (see the variable terminatef),
                  % find the numbers of function evaluations for each solver,
                  % and save them in iter(KPR, : :, :, :).
        for s = 1:KS
            for i = 1:KD*KR
                if (mod(i,KD) == 0)
                    d = KD;
                else
                    d = mod(i,KD);
                end
                r = min(floor(i/KD)+1,ceil(i/KD));
                terminatef(i)=max(max(FOPT(:,(p-1)*KD*KR+d+(0:KR-1)*KD)));
                iter(pr,s,p,d,r)=find(FREC(s,i,:) <= terminatef(i),1,'first'); 
            end
        end
    end
end

for pr = 1:KPR
    for s = 1:KS
        for p = 1:KP
            for d = 1:KD
                meaniter(pr,s,p,d)=mean(iter(pr,s,p,d,:));
                miniter(pr,s,p,d)=min(iter(pr,s,p,d,:));
                maxiter(pr,s,p,d)=max(iter(pr,s,p,d,:));
                stditer(pr,s,p,d)=std(iter(pr,s,p,d,:));
                rstditer(pr,s,p,d)=stditer(pr,s,p,d)/meaniter(pr,s,p,d);
            end
        end
    end
end
fprintf('\n');

DIMENSION=zeros(KP*KD,KS);
for i = 1:KS
    DIMENSION(:,i)=DIM;
end

clear iter FREC FRECTEMP;
save PRITERDATA  KP KPR KS KD KR DIMENSION tau meaniter miniter maxiter stditer rstditer fail;

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
%    max(T./DIMENSION)
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




