function [f, g, H] = extrosnb(x)
% A extended version of the Rosenbrock function
%
%   See CUTEst
%
%   S. Gratton, L. N. Vicente, Z. Zhang, 2019
%

alpha = 4;

n = length(x);

f = (x(1)-1)^2 + alpha*sum((x(2:end)-x(1:end-1).^2).^2);

g = zeros(n,1);
g(1) = 2*(x(1)-1);

H = zeros(n,n);
H(1,1) = 2;

for i = 2:n
    g(i) = g(i) + 2*alpha*(x(i)-x(i-1)^2);
    g(i-1) = g(i-1) + 2*alpha*(x(i)-x(i-1)^2)*(-2*x(i-1));

    H(i,i) = H(i,i) + 2*alpha;
    H(i-1, i-1) = H(i-1, i-1) + alpha*(12*x(i-1)^2-4*x(i));
    H(i-1, i) = H(i-1, i) - 4*alpha*x(i-1);
    H(i, i-1) = H(i-1, i);
end

return
