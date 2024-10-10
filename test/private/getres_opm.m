function [ rcostf, rcostg, rgapsq, errq2, errqv ] = getres_opm( prob, x, r, fval, fstar )

global costf costg model

[ fbar, gx ] = feval( prob, 'objf', x );
gap = gx - r;
if ( ischar( fstar ) || isinf(fstar) )
   rgapsq = '-';
   errq2  = '-';
   errqv  = '-';
else
   rgapsq = 0.5 * gap'*gap      / max( 1, abs( fstar ) );
   errq2  = abs( fbar - fstar ) / max( 1, abs( fstar ) );
   errqv  = abs( fval - fbar )  / max( 1, abs( fstar ) );
end
switch( model )
case 'linear_convergence'
   rcostf = costf/iter;
   rcostg = costg/iter;
case 'multiple_precision'
   rcostf = costf;
   rcostg = costg;
end

return

end
