clear all; clc;
% In this example, we use PESTO for computing the worst case value of 
%           phi_N/phi_0  (phi_i is the Lyapunov defined in the monograph)
% where xN is obtained by doing N steps of the triple momentum method,
% F is a smooth (possibly m-strongly) convex function,
% and x_* is an optimal point of F.
% (Note: by homogeneity, this ratio can be computed by fixing 
%             phi_0=R
% with R>0 some arbitrary positive constant, which we set to R=1 below).

% Parameters
L = 1;          % smoothness constant
m = .1;         % strong convexity
N = 5;          % number of iterations


% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
param.L  = L;     
param.mu = m; 
F=P.DeclareFunction('SmoothStronglyConvex',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();		 % x0 is some starting point
[xs,fs] = F.OptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
% Initial condition phi_0<= 1
[g0,f0] = F.oracle(x0);
phi0 = f0-fs-1/2/L*g0^2-m/(2*(1-m/L))*(x0-xs-1/L*g0)^2 + m/(1-m/L) * (x0-xs)^2;
P.InitialCondition( phi0 <= 1); 


% (3) Algorithm
y = x0; z = x0; q = m/L;
[gy,fy] = F.oracle(y); % evaluate gy, the gradient at y.
A(1) = 0;
for i = 1:N
    y       = (1-sqrt(q))/(1+sqrt(q)) * (y-1/L*gy) + (1- (1-sqrt(q))/(1+sqrt(q))) * z;
    [gy,fy] = F.oracle(y); % evaluate gy, the gradient at y.
    z       = sqrt(q) * (y-1/m*gy)+(1-sqrt(q)) * z;
end

% (4) Set up the performance measure
phik = fy-fs-1/2/L*gy^2-m/(2*(1-m/L))*(y-xs-1/L*gy)^2 + m/(1-m/L) * (z-xs)^2;
criterion = phik;
P.PerformanceMetric(criterion);      % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
[double(criterion) (1-sqrt(q))^(2*N)]
% the first term should be smaller than the guarantee, as it corresponds to
% the worst-case value of the ratio (but for this case, they should match
% up to numerical precision)

