function test_prodSR1

rng('default');
n = 10;
p = 5;
S = rand( n, p );
A = rand(n,n);
A = A + A';
Y = A*S;
M = eye(n);
ee = zeros(n,1);
for i = 1:n
   v = rand(n,1),
   for k = 1:p
%   k = k%D
      rk = Y(:,k)-M*S(:,k);
      if ( abs ( rk'*S(:,k ) )> max( 1e-10, 1e-12 * S(1:n,k)'*S(1:n,k) ) )
         M = M + rk * rk' / (rk'*S(:,k));
         MBS = M*[S v];
      end
   end
end
A
P = S*((S'*S)\S');
errA=norm(P*A*P-M)
for i = 1:n
   ee(i) = 1;
   AR1(1:n,i) = M*ee;
   AR2(1:n,i)=prodSR1(ee, S, Y);
   ee(i) = 0;
end
AR1
AR2
errR2= norm(AR2-AR1)