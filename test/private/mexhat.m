%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = mexhat( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The mexican hat problem with penalty parameter 0.00001
%
%   Source:
%   A.A. Brown and M. Bartholomew-Biggs,
%      "Some effective methods for unconstrained optimization based on
%      the solution of ordinary differential equations",
%      Technical Report 178, Numerical Optimization Centre, Hatfield
%      Polytechnic, (Hatfield, UK), 1987.
%
%   Ph. Toint 28 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 0.86; 0.72 ];               % x0
   varargout{2} =  0.00194568050;               % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-2-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x    = varargin{1};
   varargout = opm_eval_cpsf( 'mexhat', 'elobjf', x, { [1:2] }, nargout );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x  = varargin{2};
   kappa = 1e5;
   r1 = 2*(x(1)-1)^2;
   r2 = sqrt(kappa)*(-0.02 +kappa*(x(2)-x(1)^2) + (x(2)-1)^2);   %NEW
   varargout{1} = r1^2 + r2^2;
   if( nargout > 1 )
     J1 = [4*(x(1)-1);0];
     J2 = sqrt(kappa)*[ -2*kappa*x(1); kappa+2*(x(2)-1) ];
     varargout{2} = 2 * ( J1*r1 + J2*r2 );
     if ( nargout > 2 )
        H1 = [ 4, 0; 0, 0 ];
	H2 = sqrt(kappa)*[ -2*kappa, 0; 0, 2 ];
	varargout{3}= 2*( J1*J1.' + r1*H1 + J2*J2.' + r2*H2 );
     end
   end
end

return

end