function [f,g,H]=engval1(x)
% ENGVAL1 function
%   
%   See
%   [1] Ph.L. Toint, "Test problems for partially separable optimization
%   and results for the routine PSPMIN",  Report 83/4, Department of
%   Mathematics, FUNDP (Namur, B), 1983
%   [2] A.R. Buckley, "Test functions for unconstrained minimization",
%   TR 1989CS-3, Mathematics, statistics and computing centre, Dalhousie
%   University, Halifax (CDN), 1989
%
%   S. Gratton, L. N. Vicent, Z. Zhang, 2019

n=length(x);

f = (sum((x(1:end-1).^2+x(2:end).^2).^2-4*x(1:end-1)+3));

g = ([4*(x(1:end-1).^2+x(2:end).^2).*x(1:end-1)-4;0]...
        +[0; 4*(x(1:end-1).^2+x(2:end).^2).*x(2:end)  ]) ;

H = (spdiags( ...
[[8*x(1:end-1).*x(2:end);0] ,   [12*x(1:end-1).^2+4*x(2:end).^2;0]+[0;12*x(2:end).^2+4*x(1:end-1).^2]    , [0;8*x(1:end-1).*x(2:end)]], ...
-1:1, n, n));
