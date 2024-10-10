%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = cosine( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Another function with nontrivial groups and
%   repetitious elements.
%   NB: scaled version of COSINE.  The original version uses scal = 12.
%
%   Source:
%      N. Gould, private communication, 1997
%
%   If the dimension is unspecified, the default n = 10 is chosen.
%
%   Ph. Toint and S. Gratton, 24 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 1 )
         disp( [ ' ERROR in cosine: n = ', int2str(n),' should be > 1!' ] )
      end
   else
      n = 10;
   end
   scal  = 1;                                   % the unscaled version
   varargout{1} = ones( n, 1);                  % x0
   for i = 1:n
      pip  = exp( scal* i /(n-1) );
      varargout{1}(i) = 1 / pip;
   end
   varargout{2} = 1-n;                          % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'QUR2-AY-V-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x   = varargin{1};
   n   = length( x );
   for iel = 1:n-1
      eldom{ iel } = [ iel iel+1 ];
   end
   varargout = opm_eval_cpsf( 'cosine', 'elobjf', x, eldom, nargout, n );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel  = varargin{1};
   x    = varargin{2};
   n    = varargin{3};
   scal = 6;
   pi0  = exp( scal*( iel )/(n-1) );
   pip1 = exp( scal*(iel+1)/(n-1) );
   fiel = cos( pi0^2*x(1)^2 - 0.5*pip1*x(2) );
   varargout{1} = fiel;
   if ( nargout > 1 )
      aux  = sin( pi0^2*x(1)^2 - 0.5*pip1*x(2) );
      varargout{2} = - aux*[2*pi0^2*x(1);-0.5*pip1];
      if ( nargout > 2 )
         varargout{3} = [ -fiel*(2*pi0^2*x(1))^2-2*pi0^2*aux , 2*pi0^2*x(1)*(0.5*pip1)*fiel;
                          2*pi0^2*x(1)*(0.5*pip1)*fiel       , -(0.5*pip1)^2*fiel          ];
      end
   end

return

end