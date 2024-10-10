function  varargout = iobjf( prob, x, varargin );
%function  [ fx, [ gx, ] eff_omega ] = iobjf( prob, x, omega_f, omega_g );

%  NOTE: eff_omega is eff_omega_f when only f is required, and is eff_omega_g otherwise.
%        Ok, because, for now, eff_omega is not used.

global costf costg model

omega_f = varargin{1};

%omega_f = 1e-8;  omega_g = 1e-8;   %  for single precision runs
%omega_f = 1e-4;  omega_g = 1e-4;   %  for half precision runs

if ( nargout == 2 )
   tfx = feval( prob,'objf', x );
   if ( omega_f )
      varargout{1} = tfx + omega_f * ( 2 * rand(1,1) - 1 );
   else
      varargout{1} = tfx;
   end
elseif ( nargout == 3 )
   omega_g = varargin{2};
   [ tfx, tgx ] = feval( prob, 'objf', x );
   if ( omega_f )
      varargout{1} = tfx + omega_f * ( 2 * rand(1,1) - 1 );
   else
      varargout{1} = tfx;
   end
   if ( omega_g )
      sx    = size( x );
      onesx = ones( sx );
      done  = 0;
      omeg  = omega_g;
      ntgx  = norm( tgx );
      if ( ntgx > 1e-14 )
         while ( ~done  ) 
            gx  = tgx + omeg * ( 2 * rand( sx ) - onesx ) * ntgx;
            rerrg = norm( gx - tgx ) / norm( gx);
            done = rerrg <= omega_g+1e-15;
            omeg = 0.5 * omeg;
         end
      else
         gx = omega_g * ( 2 * rand( sx ) - onesx );
      end
      varargout{2} = gx;
   else
      varargout{2} = tgx;
   end
end

%   Compute the cost.

switch( model )

case 'linear_convergence'
   costf = costf + omega_f;
   eff_omega  = omega_f;
case 'multiple_precision'
   half_precision_gain = 4;
   if ( omega_f <= 1e-8 )
      costf     = costf + 1;
      eff_omega = 1e-16;
   elseif ( omega_f < 1e-4 )
      costf     = costf + 1 / half_precision_gain;
      eff_omega = 1e-8;
   else
      costf     = costf + 1 /  half_precision_gain^2;
      eff_omega = 1e-4;
   end
end
if ( nargout == 3 )
   switch( model )
   case 'linear_convergence'
      costg = costg + omega_g;
      eff_omega = omega_g;
   case 'multiple_precision'
      half_precision_gain = 4;
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
end
varargout{nargout} = eff_omega;

return

end
