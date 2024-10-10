function testlarge (n, signif)

plist = {'ARWHEAD', 'BRYBND',  'CHROSEN', 'cragglvy', 'DIXMAANE', 'ENGVAL1', 'EG2', 'LIARWHD', 'NONDIA', 'POWER', 'SPARSQUR', 'WOODS'};
for ip = 1 : length(plist)
    prob = plist{ip};
    %options=[]; options.maxfun = 1; options.signif=signif; testalg('newuoas', prob, n, options);
    options=[]; options.signif=signif; testalg('newuoas', prob, n, options);
    options=[]; options.signif=signif; testalg('fminunc', prob, n, options);
end
