function [p, u, v] = evaluate_SLP(xBoundary, yBoundary, xOmega, yOmega, sigma, jac)
%EVALUATE_SLP This will evaluate the single layer potential for 2D laplace.
%INPUTS
%	xBoundary, vector N, x coordinates of the boundary.
%	yBoundary, vector N, y coordinates of the boundary.
%	xOmega, real, The x coordinate to evaluate at.
%	yOmega, real, The y coordinate to evaluate at.
%	sigma, vector N, Weight function of each boundary node.
%	jac, vector N, ds/dt, arc-length scaled by the arc's angle.
%OUTPUTS
%	p, real, The pressure at the given location.
	N = length(sigma);
	p = (2*pi/N)*sum(log(sqrt((xOmega - xBoundary).^2 + (yOmega - yBoundary).^2)).*sigma.*jac)/(2*pi);
    
    u = -(2*pi/N)*sum(1./((xOmega - xBoundary).^2 + (yOmega - yBoundary).^2).*(xOmega - xBoundary).*sigma.*jac)/(2*pi);
    v = -(2*pi/N)*sum(1./((xOmega - xBoundary).^2 + (yOmega - yBoundary).^2).*(yOmega - yBoundary).*sigma.*jac)/(2*pi);