function Bp = prodSR1( p, S, Y )

%   Compute the product of the vector p with a limited memory SR1 approximation
%   of the Hessian defined by the secant pairs contained in S (steps) and Y
%   (corresponding gradient differences). The SR1 formula is given by
%
%               (y-Bs)*(y-Bs)'
%   B+  =  B +  --------------
%                 (y-Bs)'*s
%
%   Programming: Ph. Toint 1 V 2018
%

[ n, m ] = size( S );
if ( m )
   Bp = p;
   BS = S;
%   Bp = ( ( Y(1:n,1)'*(Y(1:n,1) ) / ( Y(1:n,1)'*S(1:n,1) ) ) * p;
   for k = 1:m
%       k = k%D
       rk  = Y(1:n,k) - BS(1:n,k);
       den = rk' * S(1:n,k);
       if ( abs ( den  ) > max( 1e-10, 1e-12 * S(1:n,k)'*S(1:n,k) ) )
          for j = k:m
	     BS(1:n,j) = BS(1:n,j) + rk * ( rk' * S(1:n,j) ) / den;
	  end
          Bp = Bp + rk * ( rk' * p ) / den;
       end
%       BSp = [ BS Bp]%D
   end
else
   Bp = p;
end