function [f, g, H] = minsurf(x, freq)
%A function based on the minimal surface problem defined on [0, 1] with
%   boundary conditions:
%   x(t,s) = sin(4pi*s) + 0.1*sin(120pi*s), t = 0 or 1, 0 <= s <= 1
%   x(t,s) = 0, s = 0 or 1, 0 <= t <= 1 
%   
%   Note that the length(x) should be the square of a nonzero integer
% 
%   xstart = ones(128^2,1);
%
%   S. Gratton, L. N. Vicent, Z. Zhang, 2019

% freq = 'l': Low frequency dominates in the boundary conditions
% freq = 'h': High frequency dominates in the boundary conditions
% freq = 'm': Low/high frequencies have similar magnitudes

if nargin == 1
    freq = 'l';
end

f = fn(x, freq);
g = grad(x, freq);
H = hessian(x, freq);	

function [c]=fn(MV, freq)
% [c]=cost(MV);
%  Computation of the cost of the problem
%
if size(MV,2)==1
  N=sqrt(size(MV,1));
  MV=reshape(MV,N,N);
else
  N=size(MV,1);
end
np1 = N + 1;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADD THE BOUNDARY CONDITIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
h=1/(N+1);
MM=complete(MV, freq);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STRUCTURE OF THE FUNCTION TO BE INTEGRATED 
%  THE CONTRIBUTION OF EACH PAIR OF TRIANGLE TO THE COST, 
%  GRADIENT AND HESSIAN ARE COMPUTED USING THE FUNCTION 
%  LOCAL_QUANTITIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
c=0;
% 
for i=1:np1
   for j=1:np1
      x=MM(i,j);
      y=MM(i+1,j);
      z=MM(i+1,j+1);
      t=MM(i,j+1);
      c1=sqrt(1+ ((x-t)/h)^2    + ((z-t)/h)^2);
      c2=sqrt(1+ ((x-y)/h)^2    + ((z-y)/h)^2);
      lc=(c1+c2)*h^2/2;
      c=c+lc;
   end
end


function [g]=grad(MV, freq)
% [g]=grad(MV);
% Computation of the gradient of the problem
%
if size(MV,2)==1
  N=sqrt(size(MV,1));
  MV=reshape(MV,N,N);
else
  N=size(MV,1);
end
np1 = N + 1;
np2 = N + 2;
n   = np2^2;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADD THE BOUNDARY CONDITIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
h=1/(N+1);
MM=complete(MV, freq);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STRUCTURE OF THE FUNCTION TO BE INTEGRATED 
%  THE CONTRIBUTION OF EACH PAIR OF TRIANGLE TO THE COST, 
%  GRADIENT AND HESSIAN ARE COMPUTED USING THE FUNCTION 
%  LOCAL_QUANTITIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
g=zeros(n,1);
% 
for i=1:np1
   for j=1:np1
      pt1  = (j-1) * np2 + i  ;
      pt2  = pt1+1;
      pt3  = pt2 +np2;
      pt4  = pt3-1;
%
      x=MM(i,j);
      y=MM(i+1,j);
      z=MM(i+1,j+1);
      t=MM(i,j+1);
%
      c1=sqrt(1+ ((t-x)/h)^2    + ((t-z)/h)^2);
      c2=sqrt(1+ ((x-y)/h)^2    + ((z-y)/h)^2);
      g(pt1) = g(pt1) + ( (x-t)/c1 + (x-y)/c2 )/2;
      g(pt2) = g(pt2) + (2*y-x-z)/2/c2;
      g(pt3) = g(pt3) + ( (z-t)/2/c1 + (z-y)/2/c2 );
      g(pt4) = g(pt4) + (2*t-x-z)/2/c1;
%
   end
end
% TRUNCATION OF THE GRADIENT TO REMOVE THE BOUNDARY POINTS;

ind=1:n;
ind=reshape(ind,np2,np2);
ind=ind(2:np1,2:np1);ind=ind(:);
g=g(ind);

function [H]=hessian(MV, freq)
% [H]=hessian(MV); 
% COMPUTATION OF THE HESSIAN OF THE PROBLEM 
% n=15^2;
% X=rand(n,1);
% dd=1e-8*rand(n,1);
% Xp=X+dd;
% c=cost(X);cp=cost(Xp);
% g=grad(X);gp=grad(Xp);
% H=hessian(X);
% [cp-c  ;g'*dd ]
% [gp-g,H*dd]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IF MV IS A VECTOR, TRANSFORM IT INTO A MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if size(MV,2)==1
  N=sqrt(size(MV,1));
  MV=reshape(MV,N,N);
else
  N=size(MV,1);
end
np1 = N + 1;
np2 = N + 2;
n   = np2^2;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADD THE BOUNDARY CONDITIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
h=1/(N+1);
MM=complete(MV, freq);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STRUCTURE OF THE FUNCTION TO BE INTEGRATED 
%  THE CONTRIBUTION OF EACH PAIR OF TRIANGLE TO THE COST, 
%  GRADIENT AND HESSIAN ARE COMPUTED USING THE FUNCTION 
%  LOCAL_QUANTITIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
   DD0 = zeros(n,1);
   DD1 = zeros(n-1,1);
   DD2 = zeros(n-np2,1);
   DD3 = zeros(n-np2-1,1);
% 
for i=1:np1
   for j=1:np1
      pt(1)=(j-1) * np2 + i  ;
      pt(2)=(j-1) * np2 + i+1;   
      pt(3)=  j   * np2 + i+1;
      pt(4)=  j   * np2 + i  ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      x=MM(i,j);
      y=MM(i+1,j);
      z=MM(i+1,j+1);
      t=MM(i,j+1);
   c1=sqrt(1+ ((t-x)/h)^2    + ((t-z)/h)^2);
   c2=sqrt(1+ ((x-y)/h)^2    + ((z-y)/h)^2);
   %
% lOCAL GRADIENT 
   g1=1/(h)^2/c1* [ x-t
               0  
               z-t
               2*t-x-z
                  ];
   g2=1/(h)^2/c2* [ x-y
               2*y-x-z
               z-y
               0   
                  ];
  xmt=x-t; zmt=z-t;dtmxmz=2*t-x-z;
  xmy=x-y; zmy=z-y;dymxmz=2*y-x-z;
% LOCAL HESSIAN
%  H1=1/c1^2/h^2* [ c1-g1(1)*(x-t),     -g1(2)*(x-t),         -g1(3)*(x-t),     -c1-g1(4)*(x-t);
%                               0,               0,                    0,              0;
%                      -g1(1)*(z-t),     -g1(2)*(z-t),        c1-g1(3)*(z-t),     -c1-g1(4)*(z-t);
%                 -c1-g1(1)*(2*t-x-z), -g1(2)*(2*t-x-z),   -c1-g1(3)*(2*t-x-z),2*c1-g1(4)*(2*t-x-z)];
%	   
%  H2=1/c2^2/h^2* [c2-g2(1)*(x-y),     -c2-g2(2)*(x-y),       -g2(3)*(x-y),     -g2(4)*(x-y); ...
%                 -c2-g2(1)*(2*y-x-z),2*c2-g2(2)*(2*y-x-z),-c2-g2(3)*(2*y-x-z), -g2(4)*(2*y-x-z) ; ...
%                    -g2(1)*(z-y),     -c2-g2(2)*(z-y),     c2-g2(3)*(z-y),     -g2(4)*(z-y); ...
%	                0,                   0,                  0,               0];
%lH=(H1+H2)*h^2/2;
%
%           DD0(pt,1)   =DD0(pt,1)   +diag(full(lH));
            DD0(pt,1) = DD0(pt,1) + [c1/c1^2/2-g1(1)*xmt/c1^2/2+c2/c2^2/2-g2(1)*xmy/c2^2/2;
             2*c2/c2^2/2-g2(2)*dymxmz/c2^2/2;
             c1/c1^2/2-g1(3)*zmt/c1^2/2+c2/c2^2/2-g2(3)*zmy/c2^2/2;
             2*c1/c1^2/2-g1(4)*dtmxmz/c1^2/2];
%           DD1(pt(1),1)=DD1(pt(1),1)+(lH(1,2));
            DD1(pt(1),1)=DD1(pt(1),1)-g1(2)*xmt/c1^2/2-c2/c2^2/2-g2(2)*xmy/c2^2/2;
%           DD2(pt(1),1)=DD2(pt(1),1)+(lH(1,4));
            DD2(pt(1),1)=DD2(pt(1),1)-c1/c1^2/2-g1(4)*xmt/c1^2/2-g2(4)*xmy/c2^2/2;
%           DD3(pt(1),1)=DD3(pt(1),1)+(lH(1,3));
            DD3(pt(1),1)=DD3(pt(1),1)-g1(3)*xmt/c1^2/2-g2(3)*xmy/c2^2/2;
            %
%           DD2(pt(2),1)=DD2(pt(2),1)+(lH(2,3));
            DD2(pt(2),1)=DD2(pt(2),1)-c2/c2^2/2-g2(3)*dymxmz/c2^2/2;
            %
%           DD1(pt(4),1)=DD1(pt(4),1)+(lH(3,4));
             DD1(pt(4),1)=DD1(pt(4),1)-c1/c1^2/2-g1(4)*zmt/c1^2/2-g2(4)*zmy/c2^2/2;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   end
end

%H = spdiags( [[DD3;zeros(N+3,1)] [zeros(N+3,1);DD3] ], [-(N+3) N+3 ], n, n );
%H = spdiags(DD0,0,H);
%H = spdiags([0;DD1], 1,H);
%H = spdiags([DD1;0],-1,H);
%H = spdiags([zeros(np2,1);DD2], np2,H);
%H = spdiags([DD2;zeros(np2,1)],-np2,H);
H = spdiags( [[DD3;zeros(N+3,1)] [zeros(N+3,1);DD3] DD0 [0;DD1] [DD1;0],...
      [zeros(np2,1);DD2] [DD2;zeros(np2,1)]],  [-(N+3) N+3 0 1 -1 np2 -np2 ], n, n );


% TRUNCATION OF THE HESSIAN TO REMOVE THE BOUNDARY POINTS;

ind=1:n;
ind=reshape(ind,np2,np2);
ind=ind(2:np1,2:np1);ind=ind(:);
H=H(ind,ind);

function M=complete(MV, freq)
%  M=complete(MV);
%  ADD THE BOUNDARY CONDITION TO THE 
%  THE UNKNOWN MV. MV IS EITHER A VECTOR
%  OR A MATRIX.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IF MV IS A VECTOR, TRANSFORM IT INTO A MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

lf = 4*pi; % Low frequency
hf = 120*pi; % High frequency
switch freq
case 'l'
    alf = 1; % Amplitude of low frequency
    ahf = 1e-1; % Amplitude of high frequency
case 'h'
    alf = 1e-1;
    ahf = 1;
otherwise 
    alf = 0.5;
    ahf = 0.5;
end

if size(MV,2)==1
  N=sqrt(size(MV,1));
  MV=reshape(MV,N,N);
else
  N=size(MV,1);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADD THE BOUNDARY CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NUMBER OF GRID POINTS INCLUDING THE BOUNDARY CONDITION
%
Ns=N+2;Nt=N+2;
%
M(2:Nt-1,2:Ns-1) =MV;
M(1,1:Ns)        = f0((0:(Ns-1))/(Ns-1), lf, alf, hf, ahf) ;
M(1:Nt,1)        = f1((0:(Nt-1))/(Nt-1), lf, alf, hf, ahf)';
M(Nt,1:Ns)       = f2((0:(Ns-1))/(Ns-1), lf, alf, hf, ahf) ;
M(1:Nt,Ns)       = f3((0:(Nt-1))/(Nt-1), lf, alf, hf, ahf)';

function r=f0(x, lf, alf, hf, ahf)
% North boundary condition
r= alf*sin(lf*x)+ahf*sin(hf*x); 
%r=x.*(1-x);

function r=f1(x, lf, alf, hf, ahf)
% East boundary condition
r= 0*x ;

function r=f2(x, lf, alf, hf, ahf)
% South boundary condition
r= alf*sin(lf*x)+ahf*sin(hf*x); 
%r=x.*(1-x);

function r=f3(x, lf, alf, hf, ahf)
% West boundary condition
r= 0*x ;
