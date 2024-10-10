%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = spmsqrt( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The tridiagonal matrix square root problem by Nocedal and Liu.
%   seen as a nonlinear least-squares problem.

%   Source:  problem 51 (p. 94) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   If the dimension is unspecified, the default n = 10 is chosen.
%
%   Ph. Toint 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      m = (n+2)/3;
      if ( n < 1 || round( m ) ~= m )
         disp( [ ' ERROR in spmsqrt: n = ', int2str(n),' but should satisfy n >= 4, n = 3m-2!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = 0.2*sin( ([1:n].^2)' );       % x0
   varargout{2} = 1.168931804483;               % fstar (n= 16)
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AY-V-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x     = varargin{1};
   n     = length( x );
   m     = (n+2)/3;

%  Define B
   
   ib    = 0;
   B     = sparse( m, m );
   b     = sin( ([1:3*m-2].^2).' );
   rows  = cell( m, 1 );
   cols  = cell( m, 1 );
   for j = 1:m
      if ( j == 1 )
         cols{ 1 } = [ 1 2 ];
	 B( [ 1,2 ], 1 )  = b(1:2);
	 ib = 2;
      elseif ( j == m )
         cols{ m } = [ n-1 n];
	 B( [ m-1,m ], m ) = b( ib+1:ib+2 );
      else
         cols{ j } = [ 3*(j-1):3*j-1 ];
	 B( [j-1,j,j+1], j ) = b(ib+1:ib+3);
	 ib   = ib + 3;
      end
   end
   for i = 1:m
      if ( i == 1 )
         rows{ 1 } = [ 1 3 ];
         is   = 2;
      elseif ( i == m-1 )
         rows{ m-1 } = [ n-5 n-3 n-1];
      elseif ( i == m )
         rows{ m } = [ n-2 n ];
      else
         rows{ i } = [ is is+2 is+4 ];
         is   = is + 3;
      end
   end
   B = B';

   %  Define A

   A  = B*B;

   %  Define one least-square residual for each entry of the nel entries of
   %  A (column-wise)
   %  iv is the integer vector such that each product of index iel is
   %       x(eldom{iel}(iv{iel}(1:l2))).'*x(eldom{iel}(iv{iel}(l2+1:l)))
   %  where l = length(iv{iel}} and l2 = l/2.  This trick is necessary in order
   %  to avoid repetitions in eldom{iel}.
   
   nel   = 5*m-6;
   IDX   = sparse(m,m);
   ii    = 0;
   for j = 1:m
       for i = max(1,j-1):min(j+1,m)
          ii       = ii + 1;
          IDX(i,j) = ii;
%x(ii)= ii;%D
       end
   end
   eldom = cell( nel, 1 );
   iv    = cell( nel, 1 );
   iel   = 0;
   for j = 1:m
       for i = max(1,j-2):min(j+2,m)
%entry = [ i, j ]
	 ik             = [ max(1,i-1):min(i+1,m) ];
	 ikk            = full( IDX( i, ik ) );
	 kj             = [ max(1,j-1):min(j+1,m) ];
	 kkj            = full( IDX( kj, j )' );
	 inter          = intersect( ik, kj );
	 [~, idx]       = ismember( inter, ik );
	 [~, jdx]       = ismember( inter, kj );
	 ijvars         = [ ikk(idx) kkj(jdx) ];
         iel            = iel + 1;
	 eldom{ iel }   = unique( ijvars, 'stable' );
%eldom{iel}
	 [ ~, iv{iel} ] = ismember( ijvars, eldom{ iel } );
%iv{iel}%D
%pause%D
      end
   end
   varargout = opm_eval_cpsf( 'spmsqrt', 'elobjf', x, eldom, nargout, A(:), iv );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel  = varargin{1};
   x    = varargin{2};
   lx   = length( x );
   iv   = varargin{4};
   liv  = length( iv{iel} );
   liv2 = liv / 2;
   ievr = iv{iel}(1:liv2);
   ievc = iv{iel}(liv2+1:liv);
%iel%D
%row=x(ievr)'%D
%col=x(ievc)'%D
%pause%D
   riel   = varargin{3}(iel) - x(ievr).'*x(ievc);
   varargout{1} = riel^2;
   if ( nargout > 1 )
      Jiel =  zeros( lx, 1 );
      for k = 1:liv2
         Jiel(ievr(k))= Jiel(ievr(k))-x(ievc(k));
	 Jiel(ievc(k))= Jiel(ievc(k))-x(ievr(k));
      end
      varargout{2} = 2*Jiel*riel;
      if ( nargout > 2 )
         Hiel = sparse( lx, lx );
         for k = 1:liv2
	    Hiel(ievc(k),ievr(k)) = Hiel(ievc(k),ievr(k))-1;
	    Hiel(ievr(k),ievc(k)) = Hiel(ievr(k),ievc(k))-1;
	 end
         varargout{3} = 2* (Jiel*Jiel.'+riel*Hiel);
     end
   end

end

return

end
