function [x, fx, exitflag, output] = newuoas(fun, x0, options)
%
% NEWUOAs v0.3
%
% NEWUOAs seeks the solution of an unconstrained optimization problems
% without using derivatives. It is designed for solving relatively LARGE
% problems. In my test, it generally outperformed NEWUOA for problems
% with more than 20 variables. Moreover, it succeeded in solving some problems
% with more than 10,000 variables, including ARWHEAD, CHROSEN, and SPARSQUR.
%
% It is NOT recommended to use NEWUOAs on problems with at most 10
% variables. In that case, one should use NEWUOA, which can be done by setting
% options.subspace = false.
%
% The algorithm NEWUOAs was described in the PhD thesis
%
% Z. Zhang, On Derivative-free Optimization Methods (in Chinese), PhD thesis, Institute of Computational Mathematics and Scientific/Engineering Computing, Chinese Academy of Sciences, Beijing, China, April 2012
%
% An English paper is being written.
%
% This version is not intended to be released. It is only for test.
%
% All rights reserved.
%
% ZHANG Zaikun, 08/08/2016
% Department of Applied Mathematics, The Hong Kong Polytechnic University
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin == 1)
    x0 = 0;
end
x0 = x0(:); % Work with column vectors.
n = length(x0);

% Default options:
rhobeg = 1;
rhoend = 1e-6;
maxfun = 100*n;
maxiter = 50;
ftarget = -Inf;
maxsubspacedim = min(3, n);
modeltype = 'quadratic';
debug = true;
chkfunval = false;
reproduce = false;
warndim = true;
subspace = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin == 3) % Adopt user-defined options.
    if (~isa(options, 'struct'))
        error('NEWUOAs:InvalidOptions', 'options should be a structure.');
    end
    if (isfield(options,'rhobeg'))
        rhobeg = options.rhobeg;
    end
    if (isfield(options,'rhoend'))
        rhoend = max(options.rhoend, eps);
    end
    if (isfield(options,'maxfun'))
        maxfun = max(1, int32(options.maxfun));
    end
    if (isfield(options,'maxiter'))
        maxiter = int32(options.maxiter);
    end
    if (isfield(options,'maxsubspacedim'))
        % Generally, the subspace contains at least 3 directions:
        % -g_k, d_{k-1}, and p_k, where g_k is an approximate gradient
        % at the current iterate, d_{k-1} is the last successful step,
        % and p_k is a "preconditioning" direction.
        % Therefore, we enforce that maxsubspacedim >= min(n, 3).
        maxsubspacedim = min(n, max(3, int32(options.maxsubspacedim)));
    end
    if (isfield(options,'ftarget'))
        ftarget = options.ftarget;
    end
    if (isfield(options,'modeltype'))
        if (ischar(options.modeltype))
            if (strcmpi(options.modeltype, 'quadratic') || strcmpi(options.modeltype, 'q'))
                modeltype = 'quadratic';
            elseif (strcmpi(options.modeltype, 'linear') || strcmpi(options.modeltype, 'l'))
                modeltype = 'linear';
            end
        end
    end
    if (isfield(options,'debug'))
        debug = options.debug;
    end
    if (isfield(options,'chkfunval'))
        chkfunval = options.chkfunval;
    end
    if (isfield(options,'reproduce'))
        reproduce = options.reproduce;
    end
    if (isfield(options,'warndim'))
        warndim = options.warndim;
    end
    if (isfield(options,'subspace'))
        subspace = options.subspace;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check inputs.
if (~isa(fun, 'function_handle') && ~isa(fun, 'char'))
    error('NEWUOAs:InvalidFun', 'fun should be a function handle or function name.');
end

if (~isnumeric(x0) || ~isvector(x0))
    error('NEWUOAs:InvalidX0', 'x0 should be a numerical vector or scalar.')
end
x0 = x0(:);

if (~isnumeric(rhobeg) || ~isnumeric(rhoend) || rhobeg <=0 || rhoend < 0 || rhobeg < rhoend)
    error('NEWUOAs:InvalidRhobegRhoend', 'rhobeg and rhoend should be real numbers satisfying\nrhobeg >= rhoend >= 0 and rhobeg > 0.');
end

if (~isnumeric(ftarget) || ftarget ~= ftarget)
    error('NEWUOAs:InvalidFtarget', 'ftarget should be a number.')
end

if (debug == true && chkfunval == true)
    warning('NEWUOAs:Chkfunval', 'NEWUOAs is in debug mode and checks whether the function\nvalues match the iterates, which will cost extra function evaluations. \nSet options.debug=false or options.chkfunval=false to disable the check.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Work with function handels instread of function names to avoid using 'feval'.
if (isa(fun, 'char'))
    fun = str2func(fun);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (subspace == false) % Use NEWUOA directly.
    options.subspace = false;
    [x, fx, newuoa_exitflag, newuoa_output] = newuoa(fun, x0, options);
    exitflag = 5;
    output.algorithm = 'NEWUOA';
    output.message = 'The optimization problem is solved by NEWUOA directly.';
    output.newuoa_exitflag = newuoa_exitflag;
    output.newuoa_message = newuoa_output.message;
    output.funcCount = newuoa_output.funcCount;
    output.fhist = newuoa_output.fhist;
    return;
elseif (n <= 10 && warndim == true) % Subspace methods is advantageous only when the problem is large enough.
    warning('NEWUOAs:LowDimension', 'For problems with not more than 10 variables, NEWUOA is likely to\noutperform NEWUOAs, especially if the problem is smooth and noise-free.\nSet options.subspace = false to use NEWUOA.\nSet options.warndim = false to turn off this warning.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Real NEWUOAs begins.
output.algorithm = 'NEWUOAs';
exitflag = NaN;
x = x0;
fx = fun(x);
nf = 1;
INVALID_FVALUE = -exp(163)*tan(82000)*(-sin(42)+cos(1729)); % We suppose that the objective function never returns this bizarre value. It is used only in debug mode. Contact me if your objective function does return this value.
fhist = zeros(1, maxfun) + INVALID_FVALUE;
fhist(1) = fx;

if (abs(fx) == Inf || fx ~= fx)
    exitflag = -1;
    warning('NEWUOAs:InvalidFunctionValue', 'The objective function returns an NaN or inifnite value at the starting point. NEWUOAs has to terminate.');
    output.message = 'The objective function returns an NaN or infnite value at the starting point.';
    output.iterations = 0;
    output.funcCount = nf;
    output.fhist = fhist(1:nf);
    return;
end
if (fx <= ftarget)
    exitflag = 1;
    output.message = 'NEWUOAs terminates because the target function value is achieved.';
    output.iterations = 0;
    maxiter = 0;
end

h = rhobeg;
if (exist('newuoa_options', 'var') == 1)
    clear newuoa_options;
end
newuoa_options.rhobeg = rhobeg;
newuoa_options.rhoend = rhoend;
normd = rhobeg;
D = [];
smalld = 0;
halt = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iter = 1:maxiter % Begin main loop of NEWUOAs.
    while(1)
        dim = min(size(D,2)+2, maxsubspacedim);
        if (dim >= 5)
            npt = 2*dim + 1;
        else
            npt = (dim+1)*(dim+2)/2;
        end
        if (nf >= maxfun - npt - 6 || nf >= maxfun - 2*n - npt -4 && n >= 5000) % Has to be improved.
            exitflag = 3;
            output.message = 'NEWUOAs terminates because the maximal number of function evaluation is (nearly) reached.';
            halt = true;
            break;
        end
        submaxfun = min(2*n, maxfun - nf - (npt+5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the subspace.
        if (debug == true && chkfunval == true)
            funx = fun(x);
            if (funx ~= fx && (fx == fx || funx == funx))
                error('NEWUOAs:InvalidFx', 'x and fx do not not match.');
            end
        end
        %[B, xnew, fnew, g, subnf, subfhist] = def_subspace(fun, x, fx, h, rhoend, submaxfun, modeltype, D);
        if (iter == 1 && n >= 16) % At what dimension should we use lq? For how many iterations? 2, 3, or even always?
            [B, xnew, fnew, g, subnf, subfhist] = def_subspace(fun, x, fx, min(h^2, h), rhoend, n+1, 'lq', D, reproduce); % Note that D=[] when iter=1.
        else
            [B, xnew, fnew, g, subnf, subfhist] = def_subspace(fun, x, fx, h, rhoend, submaxfun, 'quadratic', D, reproduce);
        end

        fhist(nf+1:nf+subnf) = subfhist;
        nf = nf + subnf;

        if (debug == true)
            if (subnf ~= size(subfhist, 2) || (fnew ~= min([fhist(1:nf),fnew]) && (fnew == fnew || min([fhist(1:nf),fnew]) == min([fhist(1:nf),fnew]))))
                error('NEWUOAs:def_subspace:InvalidFhist', 'def_subspace returns an subfhist that does not match subnf or fnew.');
            end
            if (chkfunval == true)
                funnew = fun(xnew);
                if (funnew ~= fnew && (fnew == fnew || funnew == funnew))
                    error('NEWUOAs:def_subspace:InvalidFnew', 'def_subspace returns xnew and fxnew that do not not match.');
                end
            end
        end

        if (isnan(sum(sum(B,2)+g+xnew)+fnew))
            exitflag = -2;
            output.message = 'def_subspace returns NaN values. This implies a bug in the code.';
            error('NEWUOAs:def_subspace:NaN', 'def_subspace returns NaN values.');
        end

        dim = size(B, 2);
        if ((dim > 0 && nnz(B) == 0) || (dim == 0 && norm(g) > 1e-3*h) || dim > maxsubspacedim)
            exitflag = -3;
            output.message = 'def_subspace returns a wrong dimension of the subspace. This implies a bug in the code.';
            error('NEWUOAs:def_subspace:WrongDim', 'def_subspace returns a wrong dimension of the subspace.');
        end

        if (dim == 0 && debug == true)
            warning('NEWUOAs:def_subspace:EmptyBasis', 'def_subspace returns an empty basis.');
        end

        x = xnew;
        fx = fnew;

        if (fx <= ftarget)
            exitflag = 1;
            output.message = 'NEWUOAs terminates because the target function value is achieved.';
            halt = true;
            break;
        end
        if (dim > 0 && norm(g)*sqrt(double(2*n)/double(min(2*n, submaxfun))) > rhoend)
            break;
        elseif (h <= rhoend)
            exitflag = 0;
            output.message = 'NEWUOAs terminates because the first-order optimality is approximately achieved.';
            output.approx_1stopt = norm(g)*sqrt(double(2*n)/double(min(2*n, submaxfun)))+h;
            output.approx_1stopt = round(output.approx_1stopt, 1, 'significant');
            halt = true;
            break;
        else
            h = h/2;
        end
    end

    if (halt == true)
        break;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call NEWUOA in the subspace.
    if (dim <= 5)
        newuoa_options.npt = (dim+1)*(dim+2)/2;
    else
        newuoa_options.npt = 2*dim+1;
    end
    newuoa_options.maxfun = min(500*dim, maxfun - nf);
    newuoa_options.rhoend = max([min(rhoend, 1.0/double(2^iter)), rhoend/max(n, 50), eps]);
    newuoa_options.rhobeg = max([newuoa_options.rhoend, h, normd, 0.5*newuoa_options.rhobeg]);
    newuoa_options.debug = debug;
    newuoa_options.chkfunval = chkfunval;
    [dopt, f, newuoa_exitflag, newuoa_output] = newuoa(@(d)fun(x+B*d), zeros(dim, 1), newuoa_options);

    normd = norm(dopt);
    if (normd ~= normd || f ~= f || f > fx || (normd > 0 && f == fx))
        exitflag = -4;
        output.message = 'NEWUOA outputs some unexpected value. This implies a bug in NEWUOA.';
        error('NEWUOAs:NEWUOAFailed', 'NEWUOA failed to solve the subproblem.\nNEWUOA outputs a step with norm %.4E and a function value %.4E', normd, f);
    end

    if (normd > 0)
        dx = B*dopt;
        x = x + dx;
        fx = f;
    end

    fhist(nf+1:nf+newuoa_output.funcCount) = newuoa_output.fhist;
    nf = nf + newuoa_output.funcCount;

    if (fx <= ftarget)
        exitflag = 1;
        output.message = 'NEWUOAs terminates because the target function value is achieved.';
        break;
    end

    if (normd <= 0.1*rhoend)
        smalld = smalld + 1;
    else
        smalld = 0;
    end
    if (smalld >= 3)
        exitflag = 2;
        output.message = 'NEWUOAs terminates because the stepsize is small for 3 consecutive iterations.';
        break;
    end

    h = max([h/2, rhoend/max(n, 50), sqrt(eps)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Up date D.
% Relatively big maxsubspacedim (default value is 3) are very
% helpful for some difficult problems including CHROSEN, VARDIM and
% WOODS. Adaptive setting of maxsubspacedim has to be investigated.
% It may have to change with n or during the iterations.
    if (normd <= 0.1*rhoend)
        dx = [];
        % When dx is small, drop it directly (can be seen as a restart
        % of the algorithm.
        % This works better than using dx = dlast+dx, which means to
        % retain the last successful step, except for a small
        % modification by dx.
    end

    if (maxsubspacedim == 3)
        D = dx;
    else
        %D = [-g, dx, D(:, 1:min(size(D,2), maxsubspacedim-4))];
        D = [dx, -g, D(:, 1:min(size(D,2), maxsubspacedim-4))];
    end
    % Notice that two other directions will join later in def_subspace.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % End main loop of NEWUOAs.

if (debug == true)
    if (nf ~= sum(fhist ~= INVALID_FVALUE) || (fx ~= min([fhist(1:nf),fx]) && (fx == fx || min([fhist(1:nf),fx]) == min([fhist(1:nf),fx]))))
        error('NEWUOAs:InvalidFhist', 'NEWUOAs returns an fhist that does not match nf or fx.');
    end
    if (chkfunval == true)
        funx = fun(x);
        if (funx ~= fx && (fx == fx || funx == funx))
            error('NEWUOAs:InvalidFx', 'NEWUOAs returns x and fx that do not not match.');
        end
    end
end

if (exitflag ~= exitflag)
    exitflag = 4;
    output.message = 'NEWUOAs terminates because the maximal number of iterations is reached.';
end

output.funcCount = nf;
if (~isfield(output, 'iterations'))
    output.iterations = iter;
end
output.fhist = fhist(1:nf);
return;


function [B, xnew, fnew, g, nf, fhist] = def_subspace(fun, x, fx, h, rhoend, maxfun, modeltype, D, reproduce)
% DEF_SUBSPACE v0.1
%
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


% Only the cases with (maxfun == 2*n && modeltype = 'quadratic') or
% (maxfun == n && modeltype = 'linear') are implemented (with some
% simple randomization, parfois). In these cases, interpolation coincide with
% finite difference, yet with an adaptively chosen stepsize, which
% stabilizes the process and makes it much more robust than fminunc of
% MATLAB (finite difference+quasi-Newton when no gradient is available).
% Finite difference is not a sin as long as it is done properly.
% Interpolation is no more than generalized and advanced finite difference.
%
% TO DO: non-finite-difference interpolation models and randomization.

% NOTE: The function evaluations in this subroutine are completely parallelable.

if (reproduce == true) % reproduce = true means to reproduce the result of the last run. This should be used only for testing algorithms.
    load('seed.newuoas.mat', 'seed');
else
    rng('shuffle');
    seed = int32(abs(1e8*(2*rand(1,1)-1))); save('seed.newuoas.mat', 'seed');
end
rng(seed);

n = length(x);

switch modeltype

case 'rq' % There must be some way to improve it. Learn from lq.
    fp = NaN(n+1,1);
    randD = 2*rand(n, n+1) - 1;
    randD = sign(randD)/sqrt(double(n));
    C = [randD', ones(n, 1)/sqrt(n)];

    for i = 1 : n+1
        fp(i) = fun(x+h*randD(:, i));
        if (abs(fp(i)) == Inf || fp(i) ~= fp(i))
            warning('NEWUOAs:def_subspace:InvalidFunctionValue', 'The objective function returns an NaN or infinite value at a point with norm %.4E', norm(x+h+randD(:, i)));
        end
    end
    nf = n+1;
    gH = (C\(fp - fx))/h;
    g = gH(1:n);
    g(g~=g) = 0;
    if(nnz(abs(g) == Inf) > 0)
        g(abs(g) < Inf) = 0;
        g(g == Inf) = 1;
        g(g == -Inf) = -1;
    end

    [fnew, ind] = min(fp);
    if (fnew < fx)
        xnew = x + h*randD(:, ind);
        if (nnz(D) > 0)
            D = D + h*randD(:, ind)*ones(1, size(D, 2));
        end
    else
        fnew = fx;
        xnew = x;
    end
    fhist = fp';
    X = [D, -g];
case 'lq'
    fp = NaN(n+1,1);
    for i = 1 : n+1
        if (i <= n)
            xtmp = x;
            xtmp(i) = x(i) + h;
        else
            xtmp = x - h*ones(n,1)/sqrt(double(n));
        end
        fp(i) = fun(xtmp);
        if (abs(fp(i)) == Inf || fp(i) ~= fp(i))
            warning('NEWUOAs:def_subspace:InvalidFunctionValue', 'The objective function returns an NaN or infinite value at a point with norm %.4E', norm(xtmp));
        end
    end
    nf = n+1;
    %ff = 0.5*((fp(n+1)-fx) + sum(fp(1:n)-fx)/sqrt(double(n))); % This works, but I do not understand why.??????
    %ff = ((fp(n+1)-fx) + sum(fp(1:n)-fx)/sqrt(double(n))); %This works even better, but I do not understand why. ?????? What is the optimal coefficient, 0.5, 1, or even bigger??????
    ff = 4.0*((fp(n+1)-fx) + sum(fp(1:n)-fx)/sqrt(double(n))); %This works even better, but I do not understand why. ?????? What is the optimal coefficient, 0.5, 1, or even bigger?????? It seems that 4 is a good choice. But why?
    if (ff == ff && abs(ff) < Inf)
        g = (fp(1:n)-fx-ff)/h;
    else
        g = (fp(1:n)-fx)/h;
    end
    g(g~=g) = 0;
    if(nnz(abs(g) == Inf) > 0)
        g(abs(g) < Inf) = 0;
        g(g == Inf) = 1;
        g(g == -Inf) = -1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Why these guys do not work?
    %C = eye(n,n) + ones(n,n)/n;
    %g = (C\(fp(1:n)-fx - (fp(n+1)-fx)/sqrt(n)))/h;  %This is the least Hession norm solution, which should work well! But it does not.
    %C = eye(n,n) + ones(n,n)/sqrt(double(n));
    %g = (C\(fp(1:n) - fp(n+1)))/h;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [fnew, ind] = min(fp);
    if (fnew < fx)
        if (ind <= n)
            xnew = x;
            xnew(ind) = x(ind) + h;
        else
            xnew = x - h*ones(n,1)/sqrt(double(n));
        end
        if (nnz(D) > 0)
            %D(ind, :) = D(ind, :) - h;
            %D(ind, :) = D(ind, :) + h;
            D = D + (xnew - x)*ones(1, size(D, 2));
        end
    else
        fnew = fx;
        xnew = x;
    end
    fhist = fp';
    X = [D, -g];
case 'quadratic'
    if (maxfun >= 2*n)
        fp = NaN(n, 1);
        fn = NaN(n, 1);
        for i = 1:n
            xtmp = x;
            xtmp(i) = x(i) + h;
            fp(i) = fun(xtmp);
            if (abs(fp(i)) == Inf || fp(i) ~= fp(i))
                warning('NEWUOAs:def_subspace:InvalidFunctionValue', 'The objective function returns an NaN or infinite value at a point with norm %.4E', norm(xtmp));
            end
            xtmp(i) = x(i) - h;
            fn(i) = fun(xtmp);
            if (abs(fn(i)) == Inf || fn(i) ~= fn(i))
                warning('NEWUOAs:def_subspace:InvalidFunctionValue', 'The objective function returns an NaN or infinite or value at a point with norm %.4E', norm(xtmp));
            end
        end
        nf = 2*n;
        g = (fp-fn)/(2*h);
        H = (fp+fn-2*fx)/(h^2);
        g(g~=g) = 0;
        if(nnz(abs(g) == Inf) > 0)
            g(abs(g) < Inf) = 0;
            g(g == Inf) = 1;
            g(g == -Inf) = -1;
        end
        H(H~=H) = 0;
        if(nnz(abs(H) == Inf) > 0)
            H(abs(H) < Inf) = 0;
            H(H == Inf) = 1;
            H(H == -Inf) = -1;
        end

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
        pg = NaN(size(g));
        threshold = max(eps, 1e-6*min(1, max(abs(H))));
        pg(H >= threshold) = g(H >= threshold)./H(H >= threshold);
        pg(H < threshold) = g(H < threshold).*(-H(H < threshold)/(threshold^2)+2/threshold);
    else
        %m = floor(maxfun/2)
        m = floor(double(maxfun)/2.0E0);
        V = randn(n, m); % Has to be improved, as m might be very close to n, which can be very large.
        [Q, ~] = qr(V, 0); % Has to be improved for the same reason.
        fp = NaN(m, 1);
        fn = NaN(m, 1);
        for i = 1:m
            fp(i) = fun(x + h*Q(:,i));
            if (abs(fp(i)) == Inf || fp(i) ~= fp(i))
                warning('NEWUOAs:def_subspace:InvalidFunctionValue', 'The objective function returns an NaN or infinite value at a point with norm %.4E', norm(x+h*Q(:,i)))
            end
            fn(i) = fun(x - h*Q(:,i));
            if (abs(fn(i)) == Inf || fn(i) ~= fn(i))
                warning('NEWUOAs:def_subspace:InvalidFunctionValue', 'The objective function returns an NaN or infinite value at a point with norm %.4E', norm(x-h*Q(:,i)))
            end
        end
        nf = 2*m;
        g = (fp-fn)/(2*h);
        H = (fp+fn-2*fx)/(h^2);
        g(g~=g) = 0;
        if(nnz(abs(g) == Inf) > 0)
            g(abs(g) < Inf) = 0;
            g(g == Inf) = 1;
            g(g == -Inf) = -1;
        end
        H(H~=H) = 0;
        if(nnz(abs(H) == Inf) > 0)
            H(abs(H) < Inf) = 0;
            H(H == Inf) = 1;
            H(H == -Inf) = -1;
        end

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
        pg = NaN(size(g));
        threshold = max(eps, 1e-6*min(1, max(abs(H))));
        pg(H >= threshold) = g(H >= threshold)./H(H >= threshold);
        pg(H < threshold) = g(H < threshold).*(-H(H < threshold)/(threshold^2)+2/threshold);
        g = Q*g;
        pg = Q*pg;
    end
    fhist = [fp', fn'];
    X = [D, -g, -pg]; % The order of the directions matters. It seems working better to put D before -g.
    %X = [D, -pg -g];
    %X = [-g, -pg, D];
case 'linear'
    if (maxfun >= n)
        fp = NaN(n,1);
        for i = 1:n
            xtmp = x;
            xtmp(i) = x(i) + h;
            fp(i) = fun(xtmp);
            if (abs(fp(i)) == Inf || fp(i) ~= fp(i))
                warning('NEWUOAs:def_subspace:InvalidFunctionValue', 'The objective function returns an NaN or infinite value at a point with norm %.4E', norm(xtmp));
            end
        end
        nf = n;
        g = (fp-fx)/h;
        g(g~=g) = 0;
        if(nnz(abs(g) == Inf) > 0)
            g(abs(g) < Inf) = 0;
            g(g == Inf) = 1;
            g(g == -Inf) = -1;
        end

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
        fp = NaN(maxfun, 1);
        V = randn(n, maxfun); % Has to be improved, as maxfun might be very close to n, which can be very large.
        [Q, ~] = qr(V, 0); % Has to be improved for the same reason.
        for i = 1:maxfun
            xtmp = x + h*Q(:, i);
            fp(i) = fun(xtmp);
            if (abs(fp(i)) == Inf || fp(i) ~= fp(i))
                warning('NEWUOAs:def_subspace:InvalidFunctionValue', 'The objective function returns an NaN or infinite value at a point with norm %.4E', norm(xtmp));
            end
        end
        nf = maxfun;
        g = (fp-fx)/h;
        g(g~=g) = 0;
        if(nnz(abs(g) == Inf) > 0)
            g(abs(g) < Inf) = 0;
            g(g == Inf) = 1;
            g(g == -Inf) = -1;
        end
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
    X = [D, -g];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dimtmp = size(X, 2);
dim = 0;
for i = 1:dimtmp
    d = X(:, i);
    if (norm(d) > 1e-3*h)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This will remove all the directions that are too short or with
        % NaN values.
        dim = dim + 1;
        d = d/norm(d, Inf);
        % This normalization will help prevent Inf/NaN from happening in
        % the subsequent computations.
        % Mysteriously, the Inf-norm normalization works a bit
        % better than the 2-norm normalization (try the L-1/4
        % regularization test problems).
        X(:, dim) = d;
    end
end
X = X(:, 1:dim);
dimtmp = dim;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim = 0;
B = NaN(size(X));
for i = 1:dimtmp
    d = X(:, i);
    normd = norm(d);
    if (normd > 1e-8)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dim = dim + 1;
        d = d/normd;
        B(:, dim) = d;
        X(:, i+1:dimtmp) = X(:, i+1:dimtmp) - d*(d'*X(:, i+1:dimtmp));
    end
end
B = B(:, 1:dim);

return;
