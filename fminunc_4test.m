function [xopt, fopt, nf] = fminunc_4test(fun, x0, place_holder, tol, maxfun)
options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', 'FunctionTolerance', tol, 'StepTolerance', tol, 'OptimalityTolerance', tol, 'MaxFunctionEvaluations', maxfun, 'MaxIterations', 500, 'Display', 'none');
[xopt, fopt,~, output] = fminunc(fun, x0, options);
nf = output.funcCount;
return;
