%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = heart6ls( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Dipole model of the heart (6 x 6 version).
%   This is a nonlinear least-squares.
%
%   Source:
%      J. E. Dennis, Jr., D. M. Gay, P. A. Vu,
%      "A New Nonlinear Equations Test Problem".
%      Tech. Rep. 83-16, Dept. of Math. Sci., Rice Univ., Houston, TX
%      June 1983, revised May 1985.
%
%   Ph. Toint 27 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]
%                   a  c  t  u  v  w
   varargout{1} = [ 0; 0; 1; 1; 1; 1 ];         % x0
   varargout{2} = 2;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AY-6-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x    = varargin{1};
   eldom{1} = [ 1:6 ];
   eldom{2} = [ 1:6 ];
   eldom{3} = [ 1:6 ];
   eldom{4} = [ 1:6 ];
   eldom{5} = [ 1:6 ];
   eldom{6} = [ 1:6 ];
   varargout = opm_eval_cpsf( 'heart6ls', 'elobjf', x, eldom, nargout );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel    = varargin{1};
   x      = varargin{2};
   sum_Mx =  -0.816;
   sum_My =  -0.017;
   sum_A  =  -1.826;
   sum_B  =  -0.754;
   sum_C  =  -4.839;
   sum_D  =  -3.259;
   sum_E  = -14.023;
   sum_F  =  15.467;
   switch ( iel )
   case 1
      iv1 = [ 3 1 ];
      iv2 = [ 4 1 ];
      iv3 = [ 5 2 ];
      iv4 = [ 6 2 ];
      switch ( nargout )
      case 1
         fe1 = heart6ls( '2prod', x(iv1) );
         fe2 = heart6ls( 'vpv'  , x(iv2), sum_Mx );
         fe3 = heart6ls( '2prod', x(iv3) );
         fe4 = heart6ls( 'vpv'  , x(iv4), sum_My );
	 varargout{1} = fe1+fe2-fe3-fe4-sum_A;
      case 2
         [ fe1, ge1 ] = heart6ls( '2prod', x(iv1) );
         [ fe2, ge2 ] = heart6ls( 'vpv'  , x(iv2), sum_Mx  );
         [ fe3, ge3 ] = heart6ls( '2prod', x(iv3) );
         [ fe4, ge4 ] = heart6ls( 'vpv'  , x(iv4), sum_My );
	 varargout{1} = fe1+fe2-fe3-fe4-sum_A;
	 varargout{2} = zeros(6,1);
	 varargout{2}(iv1) = varargout{2}(iv1) + ge1;
	 varargout{2}(iv2) = varargout{2}(iv2) + ge2;
	 varargout{2}(iv3) = varargout{2}(iv3) + ge3;
	 varargout{2}(iv4) = varargout{2}(iv4) + ge4;
      case 3
         [ fe1, ge1, He1 ] = heart6ls( '2prod', x(iv1) );
         [ fe2, ge2, He2 ] = heart6ls( 'vpv'  , x(iv2), sum_Mx  );
         [ fe3, ge3, He3 ] = heart6ls( '2prod', x(iv3) );
         [ fe4, ge4, He4 ] = heart6ls( 'vpv'  , x(iv4), sum_My );
	 varargout{1} = fe1+fe2-fe3-fe4-sum_A;
	 varargout{2} = zeros(6,1);
	 varargout{2}(iv1) = varargout{2}(iv1) + ge1;
	 varargout{2}(iv2) = varargout{2}(iv2) + ge2;
	 varargout{2}(iv3) = varargout{2}(iv3) - ge3;
	 varargout{2}(iv4) = varargout{2}(iv4) - ge4;
	 varargout{3} = zeros(6,6);
	 varargout{3}(iv1,iv1) = varargout{3}(iv1,iv1) + He1;
	 varargout{3}(iv2,iv2) = varargout{3}(iv2,iv2) + He2;
	 varargout{3}(iv3,iv3) = varargout{3}(iv3,iv3) - He3;
	 varargout{3}(iv4,iv4) = varargout{3}(iv4,iv4) - He4;
      end
   case 2
      iv5 = [ 5 1 ];
      iv6 = [ 6 1 ];
      iv7 = [ 3 2 ];
      iv8 = [ 4 2 ];
      switch ( nargout )
      case 1
         fe5 = heart6ls( '2prod', x(iv5) );
         fe6 = heart6ls( 'vpv'  , x(iv6), sum_Mx );
         fe7 = heart6ls( '2prod', x(iv7) );
         fe8 = heart6ls( 'vpv'  , x(iv8), sum_My );
	 varargout{1} = fe5+fe6+fe7+fe8-sum_B;
      case 2
         [ fe5, ge5 ] = heart6ls( '2prod', x(iv5) );
         [ fe6, ge6 ] = heart6ls( 'vpv'  , x(iv6), sum_Mx  );
         [ fe7, ge7 ] = heart6ls( '2prod', x(iv7) );
         [ fe8, ge8 ] = heart6ls( 'vpv'  , x(iv8), sum_My );
	 varargout{1} = fe5+fe6+fe7+fe8-sum_B;
	 varargout{2} = zeros(6,1);
	 varargout{2}(iv5) = varargout{2}(iv5) + ge5;
	 varargout{2}(iv6) = varargout{2}(iv6) + ge6;
	 varargout{2}(iv7) = varargout{2}(iv7) + ge7;
	 varargout{2}(iv8) = varargout{2}(iv8) + ge8;
      case 3
         [ fe5, ge5, He5 ] = heart6ls( '2prod', x(iv5) );
         [ fe6, ge6, He6 ] = heart6ls( 'vpv'  , x(iv6), sum_Mx  );
         [ fe7, ge7, He7 ] = heart6ls( '2prod', x(iv7) );
         [ fe8, ge8, He8 ] = heart6ls( 'vpv'  , x(iv8), sum_My );
	 varargout{1} = fe5+fe6+fe7+fe8-sum_B;
	 varargout{2} = zeros(6,1);
	 varargout{2}(iv5) = varargout{2}(iv5) + ge5;
	 varargout{2}(iv6) = varargout{2}(iv6) + ge6;
	 varargout{2}(iv7) = varargout{2}(iv7) + ge7;
	 varargout{2}(iv8) = varargout{2}(iv8) + ge8;
	 varargout{3} = zeros(6,6);
	 varargout{3}(iv5,iv5) = varargout{3}(iv5,iv5) + He5;
	 varargout{3}(iv6,iv6) = varargout{3}(iv6,iv6) + He6;
	 varargout{3}(iv7,iv7) = varargout{3}(iv7,iv7) + He7;
	 varargout{3}(iv8,iv8) = varargout{3}(iv8,iv8) + He8;
      end
   case 3
      iv9  = [ 1 3 5 ];
      iv10 = [ 2 3 5 ];
      iv11 = [ 1 4 6 ];
      iv12 = [ 2 4 6 ];
      switch ( nargout )
      case 1
         fe9  = heart6ls( 'adfsq', x(iv9)  );
         fe10 = heart6ls( '3prod', x(iv10) );
         fe11 = heart6ls( 'pdfsq', x(iv11), sum_Mx );
         fe12 = heart6ls( 'p3prd', x(iv12), sum_My );
	 varargout{1} = fe1-2*fe2+fe3-2*fe4-sum_C;
      case 2
         [ fe9 , ge9  ] = heart6ls( 'adfsq', x(iv9 ) );
         [ fe10, ge10 ] = heart6ls( '3prod', x(iv10) );
         [ fe11, ge11 ] = heart6ls( 'pdfsq', x(iv11), sum_Mx );
         [ fe12, ge12 ] = heart6ls( 'p3prd', x(iv12), sum_My );
	 varargout{1} = fe9-2*fe10+fe11-2*fe12-sum_C;
	 varargout{2} = zeros(6,1);
	 varargout{2}(iv9 ) = varargout{2}(iv9 ) + ge9 ;
	 varargout{2}(iv10) = varargout{2}(iv10) - 2*ge10;
	 varargout{2}(iv11) = varargout{2}(iv11) + ge11;
	 varargout{2}(iv12) = varargout{2}(iv12) - 2*ge12;
      case 3
         [ fe9 , ge9 , He9  ] = heart6ls( 'adfsq', x(iv9) );
         [ fe10, ge10, He10 ] = heart6ls( '3prod', x(iv10));
         [ fe11, ge11, He11 ] = heart6ls( 'pdfsq', x(iv11), sum_Mx );
         [ fe12, ge12, He12 ] = heart6ls( 'p3prd', x(iv12), sum_My );
	 varargout{1} = fe9-2*fe10+fe11-2*fe12-sum_C;
	 varargout{2} = zeros(6,1);
	 varargout{2}(iv9 ) = varargout{2}(iv9 ) + ge9;
	 varargout{2}(iv10) = varargout{2}(iv10) - 2*ge10;
	 varargout{2}(iv11) = varargout{2}(iv11) + ge11;
	 varargout{2}(iv12) = varargout{2}(iv12) - 2*ge12;
	 varargout{3} = zeros(6,6);
	 varargout{3}(iv9 ,iv9 ) = varargout{3}(iv9 ,iv9 ) + He9;
	 varargout{3}(iv10,iv10) = varargout{3}(iv10,iv10) - 2*He10;
	 varargout{3}(iv11,iv11) = varargout{3}(iv11,iv11) + He11;
	 varargout{3}(iv12,iv12) = varargout{3}(iv12,iv12) - 2*He12;
      end
   case 4
      iv13 = [ 2 3 5 ];
      iv14 = [ 1 3 5 ];
      iv15 = [ 2 4 6 ];
      iv16 = [ 1 4 5 ];
      switch ( nargout )
      case 1
         fe13 = heart6ls( 'adfsq', x(iv13) );
         fe14 = heart6ls( '3prod', x(iv14) );
         fe15 = heart6ls( 'pdfsq', x(iv15), sum_My );
         fe16 = heart6ls( 'p3prd', x(iv16), sum_Mx );
	 varargout{1} = fe13+2*fe14+fe15+2*fe16-sum_D;
      case 2
         [ fe13, ge13 ] = heart6ls( 'adfsq', x(iv13) );
         [ fe14, ge14 ] = heart6ls( '3prod', x(iv14) );
         [ fe15, ge15 ] = heart6ls( 'pdfsq', x(iv15), sum_My );
         [ fe16, ge16 ] = heart6ls( 'p3prd', x(iv16), sum_Mx );
	 varargout{1} = fe13+2*fe14+fe15+2*fe16-sum_D;
	 varargout{2} = zeros(6,1);
	 varargout{2}(iv13) = varargout{2}(iv13) + ge13;
	 varargout{2}(iv14) = varargout{2}(iv14) + 2*ge14;
	 varargout{2}(iv15) = varargout{2}(iv15) + ge15;
	 varargout{2}(iv16) = varargout{2}(iv16) + 2*ge16;
      case 3
         [ fe13, ge13, He13 ] = heart6ls( 'adfsq', x(iv13) );
         [ fe14, ge14, He14 ] = heart6ls( '3prod', x(iv14) );
         [ fe15, ge15, He15 ] = heart6ls( 'pdfsq', x(iv15), sum_My );
         [ fe16, ge16, He16 ] = heart6ls( 'p3prd', x(iv16), sum_Mx );
	 varargout{1} = fe13+2*fe14+fe15+2*fe16-sum_D;
	 varargout{2} = zeros(6,1);
	 varargout{2}(iv13) = varargout{2}(iv13) + ge13;
	 varargout{2}(iv14) = varargout{2}(iv14) + 2*ge14;
	 varargout{2}(iv15) = varargout{2}(iv15) + ge15;
	 varargout{2}(iv16) = varargout{2}(iv16) + 2*ge16;
	 varargout{3} = zeros(6,6);
	 varargout{3}(iv13,iv13) = varargout{3}(iv13,iv13) + He13;
	 varargout{3}(iv14,iv14) = varargout{3}(iv14,iv14) + 2*He14;
	 varargout{3}(iv15,iv15) = varargout{3}(iv15,iv15) + He15;
	 varargout{3}(iv16,iv16) = varargout{3}(iv16,iv16) + 2*He16;
      end
   case 5
      iv17 = [ 1 3 5 ];
      iv18 = [ 2 5 3 ];
      iv19 = [ 1 4 6 ];
      iv20 = [ 2 6 4 ];
      switch ( nargout )
      case 1
         fe17 = heart6ls( '3dprd', x(iv17) );
         fe18 = heart6ls( '3dprd', x(iv18) );
         fe19 = heart6ls( 'd3prd', x(iv19), sum_Mx );
         fe20 = heart6ls( 'd3prd', x(iv20), sum_My );
	 varargout{1} = fe17+fe18+fe19+fe20-sum_E;
      case 2
         [ fe17, ge17 ] = heart6ls( '3dprd', x(iv17) );
         [ fe18, ge18 ] = heart6ls( '3dprd', x(iv18) );
         [ fe19, ge19 ] = heart6ls( 'd3prd', x(iv19), sum_Mx );
         [ fe20, ge20 ] = heart6ls( 'd3prd', x(iv20), sum_My );
	 varargout{1} = fe17+fe18+fe19+fe20-sum_E;
	 varargout{2} = zeros(6,1);
	 varargout{2}(iv17) = varargout{2}(iv17) + ge17;
	 varargout{2}(iv18) = varargout{2}(iv18) + ge18;
	 varargout{2}(iv19) = varargout{2}(iv19) + ge19;
	 varargout{2}(iv20) = varargout{2}(iv20) + ge20;
      case 3
         [ fe17, ge17, He17 ] = heart6ls( '3dprd', x(iv17) );
         [ fe18, ge18, He18 ] = heart6ls( '3dprd', x(iv18) );
         [ fe19, ge19, He19 ] = heart6ls( 'd3prd', x(iv19), sum_Mx );
         [ fe20, ge20, He20 ] = heart6ls( 'd3prd', x(iv20), sum_My );
	 varargout{1} = fe17+fe18+fe19+fe20-sum_E;
	 varargout{2} = zeros(6,1);
	 varargout{2}(iv17) = varargout{2}(iv17) + ge17;
	 varargout{2}(iv18) = varargout{2}(iv18) + ge18;
	 varargout{2}(iv19) = varargout{2}(iv19) + ge19;
	 varargout{2}(iv20) = varargout{2}(iv20) + ge20;
	 varargout{3} = zeros(6,6);
	 varargout{3}(iv17,iv17) = varargout{3}(iv17,iv17) + He17;
	 varargout{3}(iv18,iv18) = varargout{3}(iv18,iv18) + He18;
	 varargout{3}(iv19,iv19) = varargout{3}(iv19,iv19) + He19;
	 varargout{3}(iv20,iv20) = varargout{3}(iv20,iv20) + He20;
      end
   case 6
      iv21 = [ 2 3 5 ];
      iv22 = [ 1 5 3 ];
      iv23 = [ 2 4 6 ];
      iv24 = [ 1 6 4 ];
      switch ( nargout )
      case 1
         fe21 = heart6ls( '3dprd', x(iv21) );
         fe22 = heart6ls( '3dprd', x(iv22) );
         fe23 = heart6ls( 'd3prd', x(iv23), sum_My );
         fe24 = heart6ls( 'd3prd', x(iv24), sum_Mx );
	 varargout{1} = fe21-fe22+fe23-fe24-sum_F;
      case 2
         [ fe21, ge21 ] = heart6ls( '3dprd', x(iv21) );
         [ fe22, ge22 ] = heart6ls( '3dprd', x(iv22) );
         [ fe23, ge23 ] = heart6ls( 'd3prd', x(iv23), sum_My );
         [ fe24, ge24 ] = heart6ls( 'd3prd', x(iv24), sum_Mx );
	 varargout{1} = fe21-fe22+fe23-fe24-sum_F;
	 varargout{2} = zeros(6,1);
	 varargout{2}(iv21) = varargout{2}(iv21) + ge21;
	 varargout{2}(iv22) = varargout{2}(iv22) - ge22;
	 varargout{2}(iv23) = varargout{2}(iv23) + ge23;
	 varargout{2}(iv24) = varargout{2}(iv24) - ge24;
      case 3
         [ fe21, ge21, He21 ] = heart6ls( '3dprd', x(iv21) );
         [ fe22, ge22, He22 ] = heart6ls( '3dprd', x(iv22) );
         [ fe23, ge23, He23 ] = heart6ls( 'd3prd', x(iv23), sum_My );
         [ fe24, ge24, He24 ] = heart6ls( 'd3prd', x(iv24), sum_Mx );
	 varargout{1} = fe21-fe22+fe23-fe24-sum_F;
	 varargout{2} = zeros(6,1);
	 varargout{2}(iv21) = varargout{2}(iv21) + ge21;
	 varargout{2}(iv22) = varargout{2}(iv22) - ge22;
	 varargout{2}(iv23) = varargout{2}(iv23) + ge23;
	 varargout{2}(iv24) = varargout{2}(iv24) - ge24;
	 varargout{3} = zeros(6,6);
	 varargout{3}(iv21,iv21) = varargout{3}(iv21,iv21) + He21;
	 varargout{3}(iv22,iv22) = varargout{3}(iv22,iv22) - He22;
	 varargout{3}(iv23,iv23) = varargout{3}(iv23,iv23) + He23;
	 varargout{3}(iv24,iv24) = varargout{3}(iv24,iv24) - He24;
      end
   end

case '2prod' % varargout = [ fiel, giel, Hiel ]

   x = varargin{1};
   varargout{1} = x(1)*x(2);
   if ( nargout > 1 )
      varargout{2} = [ x(2); x(1) ];
      if ( nargout > 2 )
         varargout{3} = [ 0, 1 ;
	                  1, 0 ];
      end
   end
   
case '3prod' % varargout = [ fiel, giel, Hiel ]

   x = varargin{1};
   varargout{1} = x(1)*x(2)*x(3);
   if ( nargout > 1 )
      varargout{2} = [ x(2)*x(3); x(1)*x(3); x(1)*x(2) ];
      if ( nargout > 2 )
         varargout{3} = [  0  , x(3), x(2) ;
	                  x(3),   0 , x(1) ;
			  x(2), x(1),   0  ];
      end
   end
   
case 'vpv' % varargout = [ fiel, giel, Hiel ]

   x     = varargin{1};
   alpha = varargin{2};
   diff  = alpha - x(2);
   varargout{1} = x(1)*diff;
   if ( nargout > 1 )
      varargout{2} = [ diff; -x(1) ];
      if ( nargout > 2 )
         varargout{3} = [  0  , -1 ;
	                  -1  ,  0 ];
      end
   end

case 'adfsq' % varargout = [ fiel, giel, Hiel ]

   x     = varargin{1};
   dfsq  = x(2)^2-x(3)^2;
   varargout{1} = x(1)*dfsq;
   if ( nargout > 1 )
      twox = 2*x(1);
      varargout{2} = [ dfsq; twox*x(2); -twox*x(3)];
      if ( nargout > 2 )
         varargout{3} = [   0   ,  2*x(2), -2*x(3);
	                  2*x(2),   twox ,   0    ;
	                 -2*x(3),    0   , -twox   ];
      end
   end

case 'pdfsq' % varargout = [ fiel, giel, Hiel ]

   x     = varargin{1};
   alpha = varargin{2};
   diff  = alpha - x(1);
   dfsq  = x(2)^2-x(3)^2;
   twod  = 2*diff;
   varargout{1} = diff*dfsq;
   if ( nargout > 1 )
      varargout{2} = [ -dfsq; twod*x(2); -twod*x(3) ];
      if ( nargout > 2 )
         varargout{3} = [     0     , -2*x(2), 2*x(3);
	                  -2*x(2),     twod  ,     0    ;
			   2*x(3),       0   ,   -twod  ];
      end
   end

case 'p3prd'

   x     = varargin{1};
   alpha = varargin{2};
   diff  = alpha - x(1);
   varargout{1} = diff*x(2)*x(3);
   if ( nargout > 1 )
      varargout{2} = [ -x(2)*x(3); diff*x(3); diff*x(2) ];
      if ( nargout > 2 )
         varargout{3} = [     0, -x(3), -x(2);
	                  -x(3),   0  ,  diff;
			  -x(2), diff ,   0  ];
      end
   end
   
case '3dprd'

   x     = varargin{1};
   diff  = x(2)^2-3*x(3)^2;
   varargout{1} = diff*x(1)*x(2);
   if ( nargout > 1 )
      varargout{2} = [ x(2)*diff; diff*x(1)+2*x(1)*x(2)^2; -6*x(1)*x(2)*x(3) ];
      if ( nargout > 2 )
         varargout{3} = [     0       , diff+2*x(2)^2, -6*x(2)*x(3);
	                 diff+2*x(2)^2,  6*x(1)*x(2) , -6*x(1)*x(3);
			 -6*x(2)*x(3) , -6*x(1)*x(3) , -6*x(1)*x(2) ];
      end
   end
   
case 'd3prd'

   x     = varargin{1};
   alpha = varargin{2};
   dfsq  = x(2)^2-3*x(3)^2;
   diff  = alpha-x(1);
   varargout{1} = diff*x(2)*dfsq;
   if ( nargout > 1 )
      varargout{2} = [ -x(2)*dfsq; diff*(dfsq+2*x(2)^2); -6*x(2)*x(3)*diff ];
      if ( nargout > 2 )
         varargout{3} = [      0       , -dfsq-2*x(2)^2,  6*x(2)*x(3);
	                 -dfsq-2*x(2)^2,   6*x(2)*diff , -6*x(3)*diff;
			  6*x(2)*x(3)  ,  -6*x(3)*diff , -6*x(2)*diff ];
      end
   end
   
end

return

end