%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = hairy( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A hairy problem in two variables.  The surface defined by
%   this function has a large number of relatively sharp hills between
%   which a valley leads to the minimizer.
%   This problem contains a large number of saddle points.
%
%   Source:
%   Ph. Toint, private communication, 1989.
%
%   Dedicated to Meret Oppenheim, creator of the "furry cup" (1936).
%
%   Ph. Toint 27 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ -5;-7];                     % x0
   varargout{2} = 20;                           % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'OUR2-AN-2-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x  = varargin{1};
   varargout = opm_eval_cpsf( 'hairy', 'elobjf', x, { [ 1 2 ] }, nargout );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel    = varargin{1};
   x      = varargin{2};
   dens   = 7;
   smooth = 0.01;
   dv1    = dens*x(1);
   dv2    = dens*x(2);
   s1sq   = sin( dv1 )^2;
   c2sq   = cos( dv2 )^2;
   fur    =  s1sq*c2sq;
   v      = x(1)-x(2);
   vsq    = v*v;
   arg    = smooth+vsq;
   dcup   = sqrt(arg);
   vsq2   = x(1)^2;
   arg2   = smooth+vsq2;
   cup1   = sqrt(arg2);
   varargout{1} = 30*fur + 100*dcup + 100*cup1;
   if ( nargout > 1 )
      tdv1   = dv1 + dv1;
      tdv2   = dv2 + dv2;
      stdv1  = sin( tdv1 );
      stdv2  = sin( tdv2 );
      gfur   = [ dens * stdv1 * c2sq; -dens * s1sq * stdv2 ];
      den    = 1/dcup;
      gdcup  = v*den*[1;-1];
      den2   = 1/cup1;
      gcup1  = [ x(1)*den2;0 ];
      varargout{2} = 30*gfur+100*gdcup+100*gcup1;
      if ( nargout > 2 )
         tdl2  = 2 * dens * dens;
         Hfur  = [ tdl2 * cos( tdv1) * c2sq, -dens * dens * stdv1 * stdv2;
	          -dens * dens * stdv1 * stdv2, -tdl2 * s1sq * cos( tdv2 ) ];
         Hdcup = (1-vsq/arg)*den*[ 1, -1;
	                          -1 , 1 ];
         Hcup1 = (1-vsq2/arg2)*den2*[ 1, 0;
	                              0, 0 ];
         varargout{3} = 30*Hfur+100*Hdcup+100*Hcup1;
      end
   end
end

return

end