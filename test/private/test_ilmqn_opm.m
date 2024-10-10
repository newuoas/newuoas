function test_ilmqn( model_in )
%
%  Tests for the inexact LMQN routine with OPM problems
%
%  Ph. Toint 7 VIII 2018

global costf costg model

probla = { 'argauss', 'arglina', 'arglinb', 'arglinc', 'argtrig', 'arwhead' };
probpa = [     3    ,     10   ,     10   ,    10    ,     10   ,     10    ];
problb = { 'bard', 'bdarwhd', 'beale', 'biggs6', 'box3', 'booth', 'brkmcc','brazthing', 'brownal','brownbs','brownden','broyden3d','broydenbd'};
probpb = [    3  ,     10   ,    2   ,    6    ,    3  ,     2  ,     2   ,     2     ,    10    ,    2    ,    4     ,    10     ,      10   ];
problc = { 'chebyqad', 'cliff', 'clustr', 'cosine', 'crglvy', 'cube' };
probpc = [     10    ,     2   ,     2       2   ,    10   ,     2  ];
probld = { 'dixmaana', 'dixmaanj', 'dixon', 'dqrtic' };
probpd = [     12  ,       12    ,   10   ,    10    ];
proble = { 'edensch',     'eg2'  , 'eg2s' , 'engval1', 'engval2' };
probpe = [     5    ,      10    ,   10   ,    10    ,    10     ];
problf = { 'freuroth' };
probpf = [     4     ];
problg = { 'genhumps', 'gottfr', 'gulf' };
probpg = [     2     ,    2    ,   4    ];
problh = { 'hairy', 'helix', 'hilbert', 'himln3', 'himm25', 'himm27', 'himm28', 'himm29', 'himm30', 'himm33', 'hypcir'};
probph = [     2  ,    3   ,     10   ,    10   ,    10   ,    10   ,   10    ,    10   ,    10   ,     10  ,    2 ];
probli = { 'indef', 'integreq' };
probpi = [     5   ,     2     ];
problj = { 'jensmp' };
probpj = [     2    ];
problk = { 'kowosb' };
probpk = [     4    ];
probll = { 'lminsurf'};
probpl = [     25   ];
problm = { 'mancino', 'mexhat', 'meyer3', 'morebv', 'msqrtals', 'msqrtbls' };
probpm = [     10   ,     2   ,    3    ,    12   ,     16    ,     16     ];
probln = { 'nlminsurf' };
probpn = [     25      ];
problo = { 'osbornea', 'osborneb' };
probpo = [      5    ,     11     ];
problp = { 'penalty1', 'penalty2', 'penalty3', 'powellbs', 'powellsg', 'powellsq', 'powr' };
probpp = [      10   ,     10    ,     10    ,     2     ,      4    ,      2    ,   10   ];
problr = { 'recipe', 'rosenbr' };
probpr = [     2   ,     2     ];
probls = { 'schmvett', 'scosine', 'sisser', 'spmsqrt' };
probps = [     3   ,        2   ,    2    ,     10    ];
problt = { 'tridia', 'trigger' };
probpt = [    10   ,    7     ];
problv = { 'vardim' };
probpv = [     10   ];
problw = { 'watson', 'wmsqrtals', 'wmsqrtbls', 'woods' };
probpw = [     12  ,     16     ,     16     ,    12   ];
problz = { 'zangwil2', 'zangwil3' };
probpz = [      2  ,        3     ];

problist = { probla{:}  , problb{:}, problc{:}, probld{:}, proble{:}, problf{:},  problg{:}, problh{:}, probli{:} };
probparm = [ probpa     , probpb   , probpc   , probpd   , probpe   , probpf   ,  probpg   , probph   , probpi    ];
problist = { problist{:}, problj{:}, problk{:}, probll{:}, problm{:}, probln{:},  problo{:}, problp{:}, problr{:} };
probparm = [ probparm   , probpj   , probpk   , probpl   , probpm   , probpn   ,  probpo   , probpp   , probpr    ];
problist = { problist{:}, probls{:}, problt{:}, problv{:}, problw{:}, problz{:}};
probparm = [ probparm   , probps   , probpt   , probpv   , probpw   , probpz   ];

problist = { 'osbornea' };
probparm = [    5  ];

%model = 'linear_convergence';
%model = 'multiple_precision';
model  = model_in;

rng(0);
nsample  = 20;

show_samples = 0;
use_diary    = 0;
save_files   = 0;
test_lmqn    = 1;
test_ilmqn   = 0;

switch( model)
case 'linear_convergence'
   resdir = './reslinear';
   if ( use_diary )
      diary './reslinear/tests.log'
   end
case 'multiple_precision'
   resdir = './resmulti';
   if ( use_diary )
      diary './resmulti/tests.log'
   end
otherwise
   disp( 'Wrong model!' )
end


%  Loop on increasing epsilons

% epss = [ 1e-3, 1e-5, 1e-7 ];
% epss = [ 1e-3, 1e-5 ];
%epss = [10^(-3)];
%epss = [10^(-5)];
%epss = [10^(-7)];
epss = [10^(-12)];
   
for ieps = 1:length( epss )
   
   eps = epss( ieps );
   acc = sprintf( '%.0e', eps );

   fprintf('\n  model = %17s    epsilon = %.2e', model, eps);
   if ( test_ilmqn )
      fprintf( '   nsample = %2d\n\n', nsample );
   else
      fprintf( '\n\n' );
   end
   disp(   '                       iters   costf     costg      rgap2     errsol    errqv ' )
   disp(   ' ' )

   %  Loop on the different problems

   for iprob = 1:length( problist );
      prob = problist{ iprob };
      disp( prob )
      [ x0, fstar ] = feval( prob, 'setup', probparm(iprob) );
      n  = length( x0 );

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ILMQN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %  Apply the exact LMQN

      if ( test_lmqn )
         maxit  = max( 1000, 5*n );
%maxit  = max( 10000, 5*n );%D
         costf  = 0;
         costg  = 0;
         Delta0 = 1;
         [ xs, fvals, rs, iter, omegas ] = ilmqn_opm( prob, x0, Delta0, eps, maxit, 0  );
if(0)%D
format long
opt=fvals(end)
%xopt=xs(:,end)'
format short
end%D
         [ rcostf, rcostg, rgapsq, errq2, errqv ] = getres_opm( prob, xs(:,end), rs(:,end), fvals(end), fstar );
         if ( isnumeric( errq2 ) )
            fprintf( '  %-17s:  %6d  %.1e   %.1e    %.1e   %.1e   %.1e\n', ...
                     'LMQN', iter, rcostf, rcostg, rgapsq, errq2, errqv )
	 else
            fprintf( '  %-17s:  %6d  %.1e   %.1e       ---       ---      ---\n', ...
                     'LMQN', iter, rcostf, rcostg )
         end	 
      end
      
      %  Apply the inexact LMQN

      if ( test_ilmqn )
	 iters  = 0;
	 costfs = 0;
	 costgs = 0;
	 errgaps = 0;
	 errq2s  = 0;
	 errqvs  = 0;
	 for is = 1:nsample
            maxit  = max( 1000, 10*n );
            costf  = 0;
            costg  = 0;
            Delta0 = 1;
            [ xs, fvals, rs, iter, omegas ] = ilmqn_opm( prob, x0, Delta0, eps, maxit, 1  );
            [ rcostf, rcostg, rgapsq, errq2, errqv ] = getres_opm( prob, xs(:,end), rs(:,end), fvals(end), fstar(1) );
	    if ( show_samples )
               if ( isnumeric( errq2 ) )
                  fprintf( '  [%3d]       :  %6d  %.1e   %.1e    %.1e   %.1e   %.1e\n', ...
                            is, iter, rcostf, rcostg, rgapsq, errq2, errqv )
	       else
                  fprintf( '  %-17s:  %6d  %.1e   %.1e       ---       ---      ---\n', ...
                           is, iter, rcostf, rcostg )
               end
            end	 
            iters   = iters + iter;
	    costfs  = costfs + rcostf;
	    costgs  = costgs + rcostg;
	    errgaps = errgaps + rgapsq;
	    if ( isnumeric( errq2 ) )
	       errq2s  = errq2s + errq2;
	    else
	       errq2s  = '-';
	    end
	    errqvs  = errqvs + errqv;
         end
         if ( isnumeric( errq2 ) )
            fprintf( '  %-17s:  %6d  %.1e   %.1e    %.1e   %.1e   %.1e\n', ...
                     'iLMQN', round(iters/nsample), costfs/nsample, costgs/nsample, errgaps/nsample, errq2s/nsample, errqvs/nsample )
         else
            fprintf( '  %-17s:  %6d  %.1e   %.1e       ---       ---      ---\n', ...
                     'iLMQN', round(iters/nsample), costfs/nsample, costgs/nsample )
         end	 
      end
      clear prob
   end
end

if ( use_diary )
   diary off
end
