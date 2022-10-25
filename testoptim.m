function testoptim(fun, n, randomx0, noise)

format long;
%randrun = 10;
randrun = 1;
rhoend = 1e-6;
maxfun = 15*n;

rng('shuffle');  
global seed i_r fval_history; 
seed = 1e8*(2*rand(1,1) - 1); save('seed.testoptim.mat', 'seed');
%load('seed.testoptim.mat', 'seed');


if (strcmp('ALL', fun))
    prob = textread('problems', '%s');
else
    prob = {fun};
end
sol = textread('solvers', '%s');

for i_p = 1 : length(prob)
    display('-------------------------------------------------------');
    display(strcat(prob{i_p}, ':'));
    [x0, rhobeg] = setup(prob{i_p}, n);
%    if (abs(randomx0)+abs(noise) > 0)
%        randrun = 3;
%    else 
%        randrun = 1;
%    end
        for i_s = 1 : length(sol)
            display(strcat(sol{i_s}, ':'));
            for i_r = 1 : randrun 
                rng(int32(1e8*abs(sin(seed) + sin(double(i_r)) + sin(sum(double(prob{i_p}))) + sin(1e8*randomx0))));
                xr = x0 + randomx0*(2*rand(n,1) - 1).*max(ones(n, 1), abs(x0)); 
                fval_history = [];
                [xopt, fopt, nf] = feval(sol{i_s}, @(x)calobjfun(prob{i_p}, x, noise), xr, rhobeg, rhoend, maxfun);
                fopt1 = min(fval_history)
                fopt2 = calobjfun(prob{i_p}, xopt, 0) 
                fopt3 = calobjfun(prob{i_p}, xopt, noise) 
                fopt
                nf
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
