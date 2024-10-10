function s = lmtrs( g, Y, S, Delta, maxit, tol )

m   = size( S, 2 );
if ( m )
   for j = 1:m
   RHO(j) = 1/(S(:,j)'*Y(:,j));
   end
   gamma_inv = (S(:,1)'*Y(:,1))/(S(:,1)'*S(:,1));
   [ a, b ]  = unrolling( S, Y, gamma_inv );
   s         = mss( g, Delta, tol, 1, S, Y, a, b, RHO, gamma_inv );
else
   s         = - Delta * g / norm( g );
end
