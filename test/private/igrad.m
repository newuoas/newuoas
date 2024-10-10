function  [ gx, eff_omega ] = igrad( g, x, omega_g );

global costf costg model

tgx   = g(x);
sx    = size( x );
onesx = ones( sx );
done  = 0;
omeg  = omega_g;
while ( ~done  ) 
  gx  = tgx + omeg * ( 2 * rand( sx ) - onesx ) * norm( tgx );
  rerrg = norm( gx - tgx ) / norm( gx) ;
  done = rerrg <= omega_g+1e-15;
  omeg = 0.5 * omeg;
end

%   Compute the cost.

switch( model )
case 'linear_convergence'
   costg = costg + omega_g;
   eff_omega = omega_g;
case 'multiple_precision'
   half_precision_gain = 4;  % for multiprecision arithmetic:
                             % should be between 3 (Horowitz et al) and 5 (Matsuoka)
                             % Can be justified as depending on the square of the
                             % bit length (area on chip)
   if ( omega_g <= 1e-8 )
      costg     = costg + 1;
      eff_omega = 1e-16;
   elseif ( omega_g < 1e-4 )
      costg     = costg + 1 / half_precision_gain;
      eff_omega = 1e-8;
   else
      costg     = costg + 1 /  half_precision_gain^2;
      eff_omega = 1e-4;
   end
end
      

return

end
