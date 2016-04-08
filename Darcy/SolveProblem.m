function [ curSig ] = SolveProblem( xq, yq, nx, ny, jac, curv, dirNodes, dirVs, neuVs )
%SOLVEPROBLEM This will solve the given problem.
%INPUTS
%   xq, vector N, The x coordinates of the boundary points.
%   yq, vector N, The y coordinates of the boundary points.
%   nx, vector N, The x values of the normals of the boundary points.
%   ny, vector N, The y values of the normals of the boundary points.
%   jac, vector N, The arc length scaled by the angular distance.
%   curv, vector N, The curvature at each boundary point.
%   dirNodes, integer vector Q, The indices of the points that have
%       Dirichlet conditions.
%   dirVs, vector N, The Dirichlet data for the nodes with Dirichlet
%       conditions. Ignore for Neuman nodes.
%   neuVs, vector N, The Neuman data for the nodes with Neuman conditions.
%       Ignore for Dirichlet nodes.
%OUTPUTS
%   curSig, vector N, The weight function at each boundary point.
    
    %TODO
    
    rhs = zeros(length(xq),1);
    
    A = mixed_poisson_SLP_matrix(xq,yq,nx,ny,jac,curv,dirNodes);
    
    rhs(dirNodes) = dirVs(dirNodes);
    rhs(setxor(dirNodes, 1:length(xq))) = neuVs(setxor(dirNodes, 1:length(xq)));
    
    curSig = A\rhs;    
    
   
end

