function [xopt, fopt, nf, fhist] = newuoas(f, x0, rhobeg, rhoend, maxfun)

n = length(x0);
if (n <= 3) 
     [xopt, fopt, nf, fhist] = newuoa(f, x0, rhobeg, rhoend, maxfun);
     return;
 end

maxiter = 50; 

x = x0;
fx = feval(f, x);
nf = 1;
fhist = fx;
h = rhobeg;
subrhobeg = rhobeg;
normd = rhobeg;
%modeltype = 'linear';
modeltype = 'quadratic';
D=[];
smalld = 0;
halt = 0;

for iter = 1:maxiter
    while(1)
        if (nf >= maxfun - 1) 
            display 1;
            halt = 1;
            break;
        end
        submaxfun = min(2*n, maxfun - nf);   
        [B, dim, xnew, fnew, g, subnf, subfhist] = def_subspace(f, x, fx, h, submaxfun, modeltype, D);
        x = xnew;
        fx = fnew;
        nf = nf + subnf;
        fhist = [fhist, subfhist];
        if (norm(g)*sqrt(double(2*n)/double(min(2*n, submaxfun))) >= rhoend)
            break;
        elseif (h <= rhoend)
            display 2;
            halt = 1;
            break;
        else
            h = h/2;
        end
    end

    if (halt == 1 || nf >= maxfun) 
            display 3;
        break;
    end

%    normg = norm(g)
%    xstar = ones(n,1);
%    in = norm(B*B'*(xstar-x))
%    out = norm(xstar-x-B*B'*(xstar-x))

    %submaxfun = maxfun - nf;
    submaxfun = min(500*dim, maxfun - nf);
    subrhoend = max(min(rhoend, 0.5^iter), rhoend/max(n, 50)); 
    subrhobeg = max([subrhoend, h, normd, 0.5*subrhobeg]);
    [dopt, fval, subnf, subfhist] = newuoa(@(d)f(x+B*d), zeros(dim, 1), subrhobeg, subrhoend, submaxfun); %!!!!!

%    normd = norm(dopt)
%    reduc = fx-fval
%    subnf
    
    x = x + B*dopt;
    fx = fval;
    nf = nf + subnf; 
    fhist = [fhist, subfhist];
    h = max([h/2, rhoend/max(n, 50), sqrt(eps)]);

    normd = norm(dopt);

    %m = 20;
    m = 10;
    %m = 1;
    if (size(D,2) <= 2*m-2)
        D = [D, -g, B*dopt];
    else
        D = [D(:,3:2*m-2), -g, B*dopt];
    end
    
    if (normd < 0.1*rhoend)
        smalld = smalld + 1;
    else
        smalld = 0;
    end
    if (smalld >= 3)
            display 4;
        break;
    end
end

xopt = x;
fopt = fx;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %tolerance = max(min(rhoend, 0.5^iter), rhoend/n);
    %tolerance = 1e-16;
    %fminunc_options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'MaxFunctionEvaluations', submaxfun, 'OptimalityTolerance', tolerance, 'StepTolerance', tolerance); 

    %[dopt, fval,~,fminunc_output] = fminunc(@(d)f(x+B*d), zeros(dim, 1), fminunc_options);
