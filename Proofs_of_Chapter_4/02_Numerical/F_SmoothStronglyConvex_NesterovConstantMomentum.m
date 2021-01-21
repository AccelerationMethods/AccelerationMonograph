clear all; clc;
% In this example, we use PESTO for computing the worst case value of 
%           (F(xN)-F_*)/(F(x0)-F_*+ m/2 * ||x0-x_*||^2)
% where xN is obtained by doing N steps of an accelerated gradient method
% with constant momentum, and F is a smooth m-strongly convex function,
% and x_* is an optimal point of F.
% (Note: by homogeneity, this ratio can be computed by fixing 
%           F(x0)-F_*+ m/2 * ||x0-x_*||^2=R
% with R>0 some arbitrary positive constant, which we set to R=1 below).

% Parameters
L = 1;          % smoothness constant
m = 0.01;       % strong convexity
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

% Initial condition F(x0)-F_*+ m/2 * ||x0-x_*||^2<= 1
[g0,f0] = F.oracle(x0);
P.InitialCondition( f0-fs+m/2*(x0-xs)^2 <= 1); 

% (3) Algorithm
x = x0; z = x0; q = m/L; A(1) = 1;

for i = 1:N
    A(i+1)  = A(i)/(1-sqrt(q));
    tau     = sqrt(q)/(1-sqrt(q));
    delta   = 1/sqrt(q);
    y       = x + tau * (z-x);
    gy      = F.gradient(y); % evaluate the gradient at y.
    x       = y -1/L * gy;
    z       = (1-q*delta)*z + q*delta*y - delta/L * gy;
end

% (4) Set up the performance measure
fN = F.value(x);
criterion = fN - fs;
P.PerformanceMetric(criterion);      % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
[double(criterion) 1/A(end) (1-sqrt(q))^N]
% the first term should be smaller than the guarantee, as it corresponds to
% the worst-case value of the ratio. The third term is equal to the second
% one.

