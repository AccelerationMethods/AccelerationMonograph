clear all; clc;
% In this example, we use PESTO for computing the worst case value of 
%           (F(xN)-F_*)/||x0-x_*||^2
% where xN is obtained by doing N steps of a gradient method,
% and F is a smooth convex function, x_* is an optimal point of F.
% (Note: by homogeneity, this ratio can be computed by fixing 
%           ||x0-x_*||=R
% with R>0 some arbitrary positive constant, which we set to R=1 below).

% Parameters
L = 1;          % smoothness constant
h = 1/L;		% step size
N = 5;          % number of iterations


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
x = x0;
A = 0;
for i = 1:N
    A = A + 1;
    x = x-h*F.gradient(x);
end

% (4) Set up the performance measure
fN = F.value(x);
criterion = fN - fs;
P.PerformanceMetric(criterion);      % Worst-case evaluated as F(x)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
[double(criterion) L/2/A L/2/N]
% the first term should be smaller than the guarantee, as it corresponds to
% the worst-case value of the ratio.
% (more precisely, it should be L/2/(2*N+1); see PESTO's documentation for 
% precise references)

