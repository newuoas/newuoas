%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = brownbs( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Brown badly scaled problem in 2 variables.  This is a nonlinear
%   least-squares problems with 3 groups.
%
%   Source: Problem 4 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Alsso problem 25 in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 17 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 1; 1 ];                     % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-2-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x  = varargin{1};
   eldom{ 1 } = [ 1 ];
   eldom{ 2 } = [ 2 ];
   eldom{ 3 } = [ 1 2 ];
   varargout = opm_eval_cpsf( 'brownbs', 'elobjf', x, eldom, nargout );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel  = varargin{1};
   x    = varargin{2};
   switch( iel )
   case 1
      riel = x(1) - 10^6;
      varargout{1} = riel^2;
      if ( nargout > 1 )
         varargout{2} = 2 * riel;
	 if ( nargout > 2 )
	    varargout{3} =  2;
	 end
      end

   case 2
      riel = x(1)-2*10^(-6);
      varargout{1} = riel^2;
      if ( nargout > 1 )
         varargout{2} = 2 * riel;
	 if ( nargout > 2 )
	    varargout{3} = 2;
	 end
      end
      
   case 3
      riel = x(1)*x(2) - 2;
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ x(2); x(1) ];
         varargout{2} = 2 * Jiel * riel;
         if ( nargout > 2 )
            varargout{3} = 2 * ( Jiel*Jiel.' + riel*[0,1;1,0] );
	 end
      end
   end
end

return

end