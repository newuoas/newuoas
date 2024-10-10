function [x, f, exitflag, output] = fminunc_4test(fun, x0, options)

if (~exist('fminunc'))
    warning('fminunc seems unavailable on your computer. Forget about testing it.');
    x = x0;
    f = feval(x0);
else
    tol = options.tol;
    maxfun = options.maxfun;

                               
    % When the function value is not accurate, we have to tune the
    % finite difference step size for fminunc. Otherwise it will not
    % work at all. Therefore, this solver is noise-aware. The other
    % solvers, including NEWUOA and NEWUOAs, are not, though they still
    % outperform fminunc on noisy problems. 
    h = sqrt(eps);
    if(isfield(options, 'noise'))
        noise = options.noise;
        if (isnumeric(noise))
            if (abs(noise) > 0)
                h = sqrt(abs(noise));
            end
        elseif (isstruct(noise))
            if (abs(noise.level) > 0)
                h = sqrt(abs(noise.level));
            end
        end
    end
    
    if (isfield(options, 'prec') && (strcmpi(options.prec, 'single') || strcmpi(options.prec, 's')))
        h = max(sqrt(sqrt(eps)), h);
    end
    
    if (isfield(options, 'signif'))
        h = max(h, sqrt(5*1.0E-1^options.signif));
    end
    
    fminunc_options = optimoptions(@fminunc, 'FiniteDifferenceStepSize', h, 'Algorithm', 'quasi-newton', 'FunctionTolerance', tol, 'StepTolerance', tol, 'OptimalityTolerance', tol, 'MaxFunctionEvaluations', maxfun, 'MaxIterations', 500, 'Display', 'none');
    
    [x, f, exitflag, output] = fminunc(fun, x0, fminunc_options);

end

return;

