function opm_describe( probname, oneline, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Prints a summary of the characteristics of the optimization probname.
%
%   INPUT:
%
%   probname: the name of the considered optimization problem
%   oneline : 1 if a one line summary (excluding header) must be produced
%   varargin: possible parameters.  If specified, varargin{1} is the problem dimension.
%
%   Programming : Ph. Toint, July 2018.
%   This version: 20 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Get problem data.

if ( nargin > 2 )
   [ x0, fstar, xtype, xlower, xupper, clower, cupper, class ] = feval( probname, 'setup', varargin{:} );
else
   [ x0, fstar, xtype, xlower, xupper, clower, cupper, class ] = feval( probname, 'setup' );
end

%  Perform its analysis.

n  = length( x0 );
if ( length( xtype ) )
   nc = length( find( xtype == 'c' ) );
   ni = length( find( xtype == 'i' ) );
else
   nc = n;
   ni = 0;
end
if ( length( xlower ) )
   xupp  = find( xupper <  Inf );
   xlow  = find( xlower > -Inf );
   xfree = setdiff( [1:n], union( xlow, xupp ) );
   xboth = intersect( xlow, xupp );
   xfix  = find( abs( xupper -  xupper ) < eps );
   xboth = setdiff( xboth, xfix );
   xlow  = setdiff( xlow, union( xboth, xfix ) );
   xupp  = setdiff( xupp, union( xboth, xfix ) );
   nfree = length( xfree );
   nlow  = length( xlow );
   nupp  = length( xupp );
   nboth = length( xboth );
   nfix  = length( xfix );
else
   nfree = n;
   nlow  = 0;
   nupp  = 0;
   nboth = 0;
   nfix  = 0;
end
m  = length( clower );
me = length( find( abs( cupper - clower ) < eps ) );
mi = m - me;

%  Test the nature and derivatives of the objective function and constraints.

if ( nargin > 2 )
   [ err_gf, err_Hf, err_gc, err_Hc, objflin, conslin ] = opm_test( probname, 0, varargin{:} );
else
   [ err_gf, err_Hf, err_gc, err_Hc, objflin, conslin ] = opm_test( probname, 0 );
end

ep = 1e-12;
if ( err_gf <= 1.e-8 && err_Hf <= 1e-8 )
   if ( objflin )
      objftype = 'linear';
   else
      objftype = 'nonlin';
   end
else
   objftype = 'wrong!';
end

if ( m )
   if ( err_gc <= 1.e-8 && err_Hc <= 1e-8 )
      if ( conslin )
         constype = 'linear';
      else
         constype = 'nonlin';
      end
   else
      constype = 'wrong!';
   end
else
   constype = '  -  ';
end

%  Printout the result.

if ( ~oneline )
   disp( ' ' )
   fprintf( 'Problem name        fstar         n   nc   ni nfree  nlow  nupp nboth  nfix  m   mi   me objtype constype  classif\n' )
   disp( ' ' )
end

if ( isnumeric( fstar ) )
   if ( fstar(1) > -Inf )
      fprintf( '%-15s %+.8e %4d %4d %4d %4d  %4d  %4d %4d  %4d %4d %4d %4d  %6s  %6s %-16s\n', ...
               probname, fstar(1), n, nc, ni, nfree, nlow, nupp, nboth, nfix, m, mi, me, objftype, constype, class );
   else
      fprintf( '%-15s      -Inf       %4d %4d %4d %4d  %4d  %4d %4d  %4d %4d %4d %4d  %6s  %6s %-16s\n', ...
               probname,  n, nc, ni, nfree, nlow, nupp, nboth, nfix, m, mi, me, objftype, constype, class );
   end
else
   fprintf( '%-15s     unknown     %4d %4d %4d %4d  %4d  %4d %4d  %4d %4d %4d %4d  %6s  %6s %-16s\n', ...
            probname,  n, nc, ni, nfree, nlow, nupp, nboth, nfix, m, mi, me, objftype, constype, class );
end

return

end
