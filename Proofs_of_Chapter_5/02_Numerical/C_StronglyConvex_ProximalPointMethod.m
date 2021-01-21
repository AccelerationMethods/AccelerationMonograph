clear all; clc;
% In this example, we use PESTO for computing the worst case value of 
%           (F(xN)-F_*)/||x0-x_*||^2
% where xN is obtained by doing N steps of a proximal point method,
% and F is a closed, proper, and (possibly m-strongly) convex function,
% and x_* is an optimal point of F.
% (Note: by homogeneity, this ratio can be computed by fixing 
%           ||x0-x_*||=R
% with R>0 some arbitrary positive constant, which we set to R=1 below).


% Parameters
m       = .1;           % strong convexity
N       = 5;            % number of iterations
lambda  = ones(N,1);	% step size (possibly a function of k)


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
x = x0;
for i = 1:N
    x = proximal_step(x, F, lambda(i));
end
xN = x;

% (4) Set up the performance measure
fN = F.value(xN);
criterion = fN-fs;
P.PerformanceMetric(criterion);

% (5) Solve the PEP
P.solve()

% (6) Compare with the bound m*||x_0-x_*||^2 / 2 / (prod(1+lambda*m)-1)
[double(criterion) m/2/(prod(1+lambda*m)-1) ]
% the first term should be smaller than the guarantee, as it corresponds to
% the worst-case value of the ratio.
% [note: the true worst-case behavior is about two times faster than the
%        guarantee in the monograph; that is, the constant is rougly 
%        divided by 2 in the sublinear regime (small strong convexity),
%        and the convergence rate is roughly squared when m is larger.]

