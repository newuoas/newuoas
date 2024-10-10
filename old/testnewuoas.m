function testnewuoas(fun, nbeg, nend, nstep)

rhoend = 1e-6;

if (strcmp('ALL', fun))
    prob = textread('problems', '%s');
    for i = 1 : lenghth(prob)
        for n = nbeg:nstep:nend
            maxfun = 50*n;
            [x0, rhobeg] = setup(fun, n);
            x0 = x0 + 10*sin(n*(1:n)').*max(ones(n, 1), abs(x0));
            [xopts, fopts, nfs, fhists] = newuoas(@(x)testfun(prob{i}, x), x0, rhobeg, rhoend, maxfun);
            [xopt, fopt, nf, fhist] = newuoa(@(x)testfun(prob{i}, x), x0, rhobeg, rhoend, maxfun);
        end
    end
else
    for n = nbeg:nstep:nend
        maxfun = 50*n;
        [x0, rhobeg] = setup(fun, n);
        x0 = x0 + 10*sin(n*(1:n)').*max(ones(n, 1), abs(x0));
        [xopts, fopts, nfs, fhists] = newuoas(@(x)testfun(fun, x), x0, rhobeg, rhoend, maxfun);
        fopts, nfs
        [xopt, fopt, nf, fhist] = newuoa(@(x)testfun(fun,x), x0, rhobeg, rhoend, maxfun);
        fopt, nf
    end
end

