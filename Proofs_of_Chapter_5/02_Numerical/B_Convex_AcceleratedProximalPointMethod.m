clear all; clc;
% In this example, we use PESTO for computing the worst case value of 
%           (F(xN)-F_*)/||x0-x_*||^2
% where xN is obtained from N steps of an accelerated proximal point method,
% and F is a closed, proper, and convex function, and x_* is an optimal
% point of F.
% (Note: by homogeneity, this ratio can be computed by fixing 
%           ||x0-x_*||=R
% with R>0 some arbitrary positive constant, which we set to R=1 below).

% Parameters
N       = 5 ;           % number of iterations
lambda  = ones(N,1);	% step size (possibly a function of k)
delta   = ones(N,1);    % possibly also a function of k ( delta(i)\in[0,1])


% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
F = P.DeclareFunction('Convex');      % F is the objective function

% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();        % x0 is some starting point
[xs,fs] = F.OptimalPoint();         % xs is an optimal point, and fs=F(xs)
P.InitialCondition( (x0-xs)^2 <= 1);% Initial condition ||x0-xs||^2 <= 1

% (3) Algorithm
A   = 0;
x = x0; z = x0;
for i = 1:N
    a         = (lambda(i)+sqrt(lambda(i)^2+4*A*lambda(i)))/2;
    y         = A/(A+a) * x + a/(A+a) * z;
    [x,gx,fx] = inexact_prox(y,F,lambda(i),delta(i));
    z         = z - a * gx;
    A         = A + a;
end

% (4) Set up the performance measure
criterion = fx - fs;
P.PerformanceMetric(criterion);

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
[double(criterion) 1/2/A 2/(sum(sqrt(lambda)))^2]
% the first term should be smaller than the guarantees (second and third
% terms), as it corresponds to the worst-case value of the ratio.
% The third term is an upper bound on the second one.

