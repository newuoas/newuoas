function [xopt, fopt, nf] = fminsearch_4test(fun, x0, place_holder, tol, maxfun)
options = optimset('MaxFunEvals', maxfun, 'TolFun', tol, 'TolX', tol, 'Display', 'none');  
[xopt, fopt,~, output] = fminsearch(fun, x0, options);
nf = output.funcCount;
return;
