function [f, g, H] = fletchcr(x, alpha)
%fletchcr: 
%    CUTEst describes it as "The chained Rosenbrock function as given by
%    Fletcher.", coming from Fletcher, "An optimal positive definite update
%    for sparse Hessian matrices", SIAM J. Optimization, 1995. In (6.5) of
%    this paper, it is defined as 
%    f(x) = \sum_{i=1}^{n-1} 100(x_{i+1} - x_i^2)^2 + (1-x_i)^2.
%    Due to the very large coefficient 100, this problem is too difficult
%    (Toint once mentioned it in private communication). 
%    Powell's version (named CHROSEN, from "The NEWUOA software for
%    unconstrained optimization without derivatives", in: Di Pillo G.,
%    Roma M. (eds) Large-Scale Nonlinear Optimization. Nonconvex Optimization
%    and Its Applications, vol 83. Springer, Boston, MA, 2006) is
%    moderately difficult with virtually the same structure. In (8.10)
%    of the paper, CHROSEN is defined as 
%    f(x) = \sum_{i=1}^{n-1} 4(x_i-x_{i+1}^2)^2 + (1-x_{i+1})^2.
%    The only essential difference is the coefficient 4 versus 100.     
%    Toint (1978) "Some Numerical Results Using a Sparse Matrix Updating
%    Formula in Unconstrained Optimization" defined a similar problem
%    (CR) in Section 3.3, namely
%    f(x) = \sum_{i=2}^n 4alpha_i*(x_{i-1}-x_i^2)^2 + (1-x_i)^2,
%    where alpha_i comes from a table of (random?) numbers ranging
%    between 0 and 3. 
%    There is yet another possible variant
%    f(x) = \sum_{i=1}^{n-1} alpha*(1+sin(i))^2*(x_i-x_{i+1}^2)^2 + (1-x_{i+1})^2, 
%    with, for example, alpha = 1. 
%
%   See 
%   [1] Fletcher, "An optimal positive definite update for sparse Hessian matrices", SIAM J. Optimization, 1995. 
%
%   L. N. Vicente, S. Gratton, Z. Zhang, 2019

n=length(x);

if nargin < 2 
    alpha = 100;
end

f=0;
g=zeros(n,1);
H=zeros(n,n);

for i=1:n-1
  f = f + (x(i)-1)^2+alpha*(x(i)^2-x(i+1))^2; 
%
  g(i)   = g(i) + 2*(x(i)-1)+alpha*2*(x(i)^2-x(i+1))*2*x(i);
  g(i+1) = g(i+1) - alpha*2*(x(i)^2-x(i+1));
%
  H(i,i)    =  H(i,i)+2+alpha*2*2*(3*x(i)^2-x(i+1));
  H(i,i+1)  =  H(i,i+1)-alpha*2*2*x(i);
  H(i+1,i)  =  H(i+1,i) -alpha*2*2*x(i);
  H(i+1,i+1)=  H(i+1,i+1)+alpha*2;
end

return
