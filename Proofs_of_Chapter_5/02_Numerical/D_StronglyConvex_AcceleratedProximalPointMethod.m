clear all; clc;
% In this example, we use PESTO for computing the worst case value of 
%           (F(xN)-F_*)/||x0-x_*||^2
% where xN is obtained from N steps of an accelerated proximal point method,
% and F is a closed, proper, and (possibly m-strongly) convex function, 
% and x_* is an optimal point of F.
% (Note: by homogeneity, this ratio can be computed by fixing 
%           ||x0-x_*||=R
% with R>0 some arbitrary positive constant, which we set to R=1 below).

% Parameters
m       = .1;               % strong convexity
N       = 5 ;               % number of iterations
lambda  = ones(N,1);        % step size (possibly a function of k)
delta   = sqrt(1+m*lambda); % also possibly a function of k
                            % with delta(i) \in [0,\sqrt{1+m*lambda(i))]

% (0) Initialize an empty PEP
P = pep();

% (1) Set up the objective function
param.mu = m;
param.L  = Inf; % Inf: no smoothness (gradient is not Lipschitz)
F = P.DeclareFunction('SmoothStronglyConvex',param);  % F is the objective function

% (2) Set up the starting point and initial condition
x0      = P.StartingPoint();        % x0 is some starting point
[xs,fs] = F.OptimalPoint();         % xs is an optimal point, and fs=F(xs)
P.InitialCondition( (x0-xs)^2 <= 1);% Initial condition ||x0-xs||^2 <= 1

% (3) Algorithm
A(1)   = 0;
x = x0; z = x0;
for i = 1:N
    A(i+1)    = A(i) + (lambda(i)+2*A(i)*lambda(i)*m+sqrt(4*A(i)^2*lambda(i)*m*(m*lambda(i)+1)+4*A(i)*lambda(i)*(lambda(i)*m+1)+lambda(i)^2))/2;
    y         = x + (A(i+1)-A(i))*(1+m*A(i))/(A(i+1)+2*m*A(i)*A(i+1)-m*A(i)^2)*(z-x);
    [x,gx,fx] = inexact_prox(y,F,lambda(i),delta(i));
    z         = z + m* (A(i+1)-A(i))/(1+m*A(i+1))*(x-z)- (A(i+1)-A(i))/(1+m*A(i+1)) * gx;
end

% (4) Set up the performance measure
criterion = fx - fs;
P.PerformanceMetric(criterion);

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
[double(criterion) 1/2/A(end) 1/2*min([4/(sum(sqrt(lambda)))^2, prod((1-sqrt(m*lambda(2:end)./(1+m*lambda(2:end)))))/lambda(1)])]
% the first term should be smaller than the guarantees (second and third
% terms), as it corresponds to the worst-case value of the ratio.
% The third term is an upper bound on the second one.

