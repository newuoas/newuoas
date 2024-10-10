function verify(solver, fun, n)

format long;
randrun = 5;
rhoend = 1e-6;

rng('shuffle');
global seed i_r fval_history;
seed = 1e8*(2*rand(1,1) - 1); save('seed.testoptim.mat', 'seed');
%load('seed.testoptim.mat', 'seed');

if (strcmp('ALL', fun))
    prob = textread('problems', '%s');
else
    prob = {fun};
end
if (strcmp('ALL', solver))
    sol = textread('solvers', '%s');
else
    sol = {solver};
end

for i_p = 1 : length(prob)
    display('-------------------------------------------------------');
    display(strcat(prob{i_p}, ':'));
    [x0, rhobeg] = setup(prob{i_p}, n);
        for i_s = 1 : length(sol)
            display(strcat(sol{i_s}, ':'));
            for i_r = 1 : randrun
                rng(int32(1e8*abs(sin(seed) + sin(double(i_r)) + sin(sum(double(prob{i_p}))))));
                xr = x0 + (2*rand(n,1) - 1).*max(ones(n, 1), abs(x0));
                for maxfun = [0, 1, n, 2*n, 2*n+1, 2*n+2, 3*n, (n+1)*(n+2)/2-1, (n+1)*(n+2)/2, (n+1)*(n+2)/2+1, n^2]
                    fval_history = [];
                    [xopt, fopt, nf] = feval(sol{i_s}, @(x)calobjfun(prob{i_p}, x, 0), xr, rhobeg, rhoend, maxfun);
                    fopt1 = calobjfun(prob{i_p}, xopt, 0);
                    fopt2 = min(fval_history);
                    if (fopt1 ~= fopt || fopt2 ~=fopt)
                        display Fail!!
                        fopt
                        fopt1
                        fopt2
                    else
                        display Pass!
                    end
                end
            end
        end
end

return;

function f = calobjfun(fun, x, noise)
global seed i_r fval_history;

f = testfun(fun, x);
fval_history = [fval_history, f];

if (isnumeric(noise))
    if (abs(noise) > 0)
        rng(int32(1e8*abs(sin(seed) + sin(double(i_r)) + sin(sum(abs(sin(1e8*x)))) + sin(sum(double(fun))) + sin(1e8*noise))));
        f = f + abs(f)*noise*(2*rand(1,1)-1);
    end
elseif (isstruct(noise))
    if (abs(noise.level) > 0)
        if (strcmpi(noise.type, 'relative') || strcmpi(noise.type, 'multiplicative') || strcmpi(noise.type, 'multiply') || strcmpi(noise.type, 'multip') || strcmpi(noise.type, '*'))
            rng(int32(1e8*abs(sin(seed) + sin(double(i_r)) + sin(sum(abs(sin(1e8*x)))) + sin(sum(double(fun))) + sin(1e8*noise.level))));
            f = f + abs(f)*noise.level*(2*rand(1,1)-1);
        elseif (strcmpi(noise.type, 'absolute') || strcmpi(noise.type, 'additive') || strcmpi(noise.type, 'add') || strcmpi(noise.type, '+'))
            rng(int32(1e8*abs(sin(seed) + sin(double(i_r)) + sin(sum(abs(sin(1e8*x)))) + sin(sum(double(fun))) + sin(1e8*noise.level))));
            f = f + noise.level*(2*rand(1,1)-1);
        end
    end
end
return;
