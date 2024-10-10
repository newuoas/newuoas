function [f, g, H, fopt, xopt] = laplace(x)
%A quadratic function based on a (perturbed) Laplace equation
%   Laplace u(x) = sin(pi x), 0 <= x <= 1
%   u(0) = u(1) = 0.
%   
%   S. Gratton, L. N Vicente, Z. Zhang, 2019

n=length(x);

h = 1/(n+1);
e = ones(n,1);
%epsilon = 5e-4; % A perturbation to Laplace
epsilon = 0; % A perturbation to Laplace
A = (1/h^2)*spdiags([-e (2+epsilon)*e -e], -1:1, n, n);

xopt = sin(10*pi*(1:n)*h)';
b   = A*xopt;
fopt = (0.5*xopt'*A*xopt-b'*xopt);

f = (0.5*x'*A*x-b'*x);
g = (A*x-b);
H=A;

