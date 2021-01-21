clear all; clc;
% In this example, we use PESTO for computing the worst case value of 
%           (F(xN)-F_*)/||x0-x_*||^2
% where xN is obtained by doing N steps of the optimized gradient method,
% and F is a smooth convex function, and x_* is an optimal point of F.
% (Note: by homogeneity, this ratio can be computed by fixing 
%           ||x0-x_*||=R
% with R>0 some arbitrary positive constant, which we set to R=1 below).

% Parameters
L = 1;          % smoothness constant
N = 5;          % number of iterations


% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
param.L  = L;     % Smoothness parameter
F = P.DeclareFunction('SmoothStronglyConvex',param); % F is the objective function

% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();		 % x0 is some starting point
[xs,fs] = F.OptimalPoint(); 		 % xs is an optimal point, and fs=F(xs)
P.InitialCondition( (x0-xs)^2 <= 1); % Initial condition ||x0-xs||^2<= 1

% (3) Algorithm
x        = cell(N+1,1);% store iterates in a cell
x{1}     = x0;
y        = x0;
theta    = cell(N+1,1);
theta{1} = 1;
for i = 1:N
    x{i+1} = gradient_step(y,F,1/param.L);
    if i<N
        theta{i+1}  = (1+sqrt(4*theta{i}^2+1))/2;
    else
        theta{i+1}  = (1+sqrt(8*theta{i}^2+1))/2;
    end
    y = x{i+1}+(theta{i}-1)/theta{i+1}*(x{i+1}-x{i})+...
        theta{i}/theta{i+1}*(x{i+1}-y);
end

% (4) Set up the performance measure
fN = F.value(y);
criterion = fN - fs;              
P.PerformanceMetric(criterion);      % Worst-case evaluated as F(y)-F(xs)

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
[double(criterion)  L/2/theta{N+1}^2  2*L/(N+2)^2]  
% the first term should be smaller than the guarantee, as it corresponds to
% the worst-case value of the ratio. The second term corresponds to the
% guarantee and the lower complexity bound, so it should match the exact
% worst-case bound provided by the first term.
% The third term is an upper bound.
