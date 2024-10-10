function [B, dim, xnew, fnew, g, nf, fhist] = def_subspace(f, x, fx, h, maxfun, modeltype, D)
% subspace is a function to generate a subspace to be used in the
% minimization of the objective function f. 
%
% Inputs:
% f: objective function handle
% x: the base point 
% fx:
% h: 
% maxfun:
% modeltype:
% D:

rng('shuffle'); 
seed = int32(abs(1e8*(2*rand(1,1) - 1))); save('seed.newuoas.mat', 'seed'); rng(seed);
%load('seed.newuoas.mat', 'seed'); rng(seed);

n = length(x);

switch modeltype
case 'quadratic'
    if (maxfun >= 2*n) 
        fp = zeros(n, 1) + NaN;
        fn = zeros(n, 1) + NaN;
        for i = 1:n
            xtmp = x;
            xtmp(i) = x(i) + h;
            fp(i) = feval(f, xtmp);
            xtmp(i) = x(i) - h;
            fn(i) = feval(f, xtmp);
        end
        nf = 2*n;
        g = (fp-fn)/(2*h);
        H = (fp+fn-2*fx)/(h^2);
        [fnew, ind] = min([fp', fn']);
        if (fnew < fx)
            if (ind <= n)
                xnew = x;
                xnew(ind) = x(ind) + h;
                g(ind) = g(ind) + h*H(ind);
                if (nnz(D) > 0) 
                    %D(ind, :) = D(ind, :) - h;
                    D(ind, :) = D(ind, :) + h;
                end
            else
                xnew = x;
                xnew(ind-n) = x(ind-n) - h;
                g(ind-n) = g(ind-n) - h*H(ind-n);
                if (nnz(D) > 0) 
                    %D(ind-n, :) = D(ind-n, :) + h;
                    D(ind-n, :) = D(ind-n, :) - h;
                end
            end
        else
            fnew = fx;
            xnew = x;
        end
        pg = g + NaN;
        threshold = 1e-6*min(1, max(abs(H)));
        ind = find(H >= threshold);
        pg(ind) = g(ind)./H(ind);
        ind = find(H < threshold);
        pg(ind) = g(ind).*(-H(ind)/(threshold^2)+2/threshold);
    else
        %m = floor(maxfun/2)
        m = floor(double(maxfun)/2.0E0);
        V = randn(n, m);
        [Q, ~] = qr(V, 0);
        fp = zeros(m, 1) + NaN;
        fn = zeros(m, 1) + NaN;
        for i = 1:m
            fp(i) = feval(f, x + h*Q(:,i));
            fn(i) = feval(f, x - h*Q(:,i));
        end
        nf = 2*m;
        g = (fp-fn)/(2*h);
        H = (fp+fn-2*fx)/(h^2);
        [fnew, ind] = min([fp', fn']);
        if (fnew < fx)
            if (ind <= m)
                xnew = x + h*Q(:, ind);
                g(ind) = g(ind) + h*H(ind);
                if (nnz(D) > 0)
                    %D = D - h*Q(:, ind)*ones(1, size(D, 2));
                    D = D + h*Q(:, ind)*ones(1, size(D, 2));
                end
            else
                xnew = x - h*Q(:, ind-m);
                g(ind-m) = g(ind-m) - h*H(ind-m);
                if (nnz(D) > 0)
                    %D = D + h*Q(:, ind-m)*ones(1, size(D, 2));
                    D = D - h*Q(:, ind-m)*ones(1, size(D, 2));
                end
            end
        else
            fnew = fx;
            xnew = x;
        end
        pg = g + NaN;
        threshold = 1e-6*min(1, max(abs(H)));
        ind = find(H >= threshold);
        pg(ind) = g(ind)./H(ind);
        ind = find(H < threshold);
        pg(ind) = g(ind).*(-H(ind)/(threshold^2)+2/threshold);
        g = Q*g;
        pg = Q*pg;
    end
    fhist = [fp', fn'];
    X = [-g, -pg, D];
case 'linear'
    if (maxfun >= n)
        fp = zeros(n,1) + NaN;
        for i = 1:n
            xtmp = x;
            xtmp(i) = x(i) + h;
            fp(i) = feval(f, xtmp);
        end
        nf = n;
        g = (fp-fx)/h;
        [fnew, ind] = min(fp);
        if (fnew < fx)
            xnew = x;
            xnew(ind) = x(ind) + h;
            if (nnz(D) > 0) 
                %D(ind, :) = D(ind, :) - h;
                D(ind, :) = D(ind, :) + h;
            end
        else
            fnew = fx;
            xnew = x;
        end
    else
        fp = zeros(maxfun, 1) + NaN;
        V = randn(n, maxfun);
        [Q, ~] = qr(V, 0);
        for i = 1:maxfun
            xtmp = x + h*Q(:, i);
            fp(i) = feval(f, xtmp);
        end
        nf = maxfun;
        g = (fp-fx)/h;
        g = Q*g;
        [fnew, ind] = min(fp);
        if (fnew < fx)
            xnew = x + h*Q(:, ind);
            if (nnz(D) > 0) 
                %D = D - h*Q(:, ind)*ones(1, size(D,2));
                D = D + h*Q(:, ind)*ones(1, size(D,2));
            end
        else
            fnew = fx;
            xnew = x;
        end
    end
    fhist = fp';
    X = [-g, D];
end

dimtmp = size(X, 2);
dim = 0;
B = X+NaN;
for i = 1:dimtmp
    d = X(:, i);
    normd = norm(d);
    if (normd > 1e-6*h)
        dim = dim + 1;
        d = d/normd;
        B(:, dim) = d;
        X(:, i+1:dimtmp) = X(:, i+1:dimtmp) - d*(d'*X(:, i+1:dimtmp)); 
        %!!! DO NOT use d*d'*X(:, i+1:dimtmp) to replace the last term;
        %!!! Otherewise, we would have a n*n >> n*timtmp dimensional
        %!!! matrix to deal with! 
    end
end
B = B(:, 1:dim);
