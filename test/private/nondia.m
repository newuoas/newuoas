function [f, g, H] = nondia(x)
% Nondiagonal variant of the Rosenbrock function
%   
%   See
%   [1] D. F. Shanno, On Variable-Metric Methods for Sparse Hessians, 1980
%
%   S. Gratton, L. N. Vicente, Z. Zhang, 2019
%

alpha = 100;

n = length(x);

f = (x(1)-1)^2 + alpha*sum((x(2:end).^2-x(1)).^2);

g = zeros(n,1);
g(1) = 2*(x(1)-1);

H = zeros(n,n);
H(1,1) = 2;

for i = 2:n
    g(1) = g(1) + 2*alpha*(x(1)-x(i)^2);
    g(i) = g(i) + 2*alpha*(x(1)-x(i)^2)*(-2*x(i));

    H(1,1) = H(1,1) + 2*alpha;
    H(i, i) = H(i, i) + alpha*(12*x(i)^2-4*x(1));
    H(1, i) = H(1, i) - 4*alpha*x(i);
    H(i, 1) = H(1, i);
end

return
