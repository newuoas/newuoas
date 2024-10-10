%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = himm25( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A simple quadratic.
%
%   Source: problem 215 (p. 61) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 20 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 0; 2 ];                     % x0
   varargout{2} = 0;                            % fstar
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-2-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x   = varargin{1};
   varargout = opm_eval_cpsf( 'himm25', 'elobjf', x, { [ 1 2 ] }, nargout);

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x  = varargin{2};
   r1 = 2 * (x(1)-5);
   r2 = (x(2) - 6 );
   varargout{1} = r1^2 + r2^2;
   if ( nargout > 1 )
      J1 = [ 2; 0 ];
      J2 = [ 0; 1 ];
      varargout{2} = 2 * ( J1*r1 + J2*r2 );
      if ( nargout > 2 )
         varargout{3} = 2 * ( J1*J1.' + J2*J2.' );
      end
   end
end

return

end