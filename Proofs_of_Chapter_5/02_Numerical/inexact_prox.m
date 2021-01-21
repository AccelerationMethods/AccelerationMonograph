function [x,gx,fx] = inexact_prox(x0,F,lambda,delta)
% This function implements an inexact proximal operator for the Performance Estimation Toolbox (PESTO)
% that is in line with the inexactness model used in the monograph.
%
% That is, the function outputs an approximate proximal solution [x,gx,fx] satisfying
%   ||x-x0+lambda*gx|| <= delta * || x - x0 ||,
% where gx is a subgradient of F at x.

% Create the output
x   = Point('Point'); 
gx  = Point('Point'); 
fx  = Point('Function value');

% Add the constraints on the input/output using PESTO functions.
% (1) associate (x,gx,fx) to F, which constrain gx to be a subgradient of F at x, and fx=F(x).
F.AddComponent(x,gx,fx); % gx is a subgradient of F at x, and fx=F(x).
% (2) add the constraint on the quality of the proximal solution.
F.AddConstraint( (x-x0+lambda*gx)^2 <= delta^2 * (x-x0)^2);

end

