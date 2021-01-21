clear all; clc;
% In this example, we use PESTO for computing the worst case value of 
%           (F(xN)-F_*)/||x0-x_*||^2
% where xN is obtained by doing N steps of an accelerated gradient method,
% and F is a smooth convex function, and x_* is an optimal point of F.
% (Note: by homogeneity, this ratio can be computed by fixing 
%           ||x0-x_*||=R
% with R>0 some arbitrary positive constant, which we set to R=1 below).

% Parameters
L = 1;          % smoothness constant
N = 2;          % number of iterations


% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
param.L  = L;     % Smoothness parameter
F=P.DeclareFunction('SmoothStronglyConvex',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();		 % x0 is some starting point
[xs,fs] = F.OptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
P.InitialCondition( (x0-xs)^2 <= 1); % Initial condition ||x0-xs||^2<= 1

% (3) Algorithm
x = x0; z = x0;
A(1) = 0;
for i = 1:N
    a       = 1/2 * (1+sqrt(4*A(i)+1));
    A(i+1)  = a + A(i);
    y       = x + (1-A(i)/A(i+1)) * (z-x);
    gy      = F.gradient(y); % evaluate the gradient at y.
    x       = y -1/L * gy;
    z       = z - (A(i+1)-A(i))/L * gy;
end

% (4) Set up the performance measure
fN = F.value(x);
criterion = fN - fs;
P.PerformanceMetric(criterion);      % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
[double(criterion) L/2/A(end) 2*L/N^2]
% the first term should be smaller than the guarantee, as it corresponds to
% the worst-case value of the ratio.

