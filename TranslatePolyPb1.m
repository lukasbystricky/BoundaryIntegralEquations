function [jac, dirVs, neuVs, dirNodes] = TranslatePolyPb1(xq, yq, markQ, nx, ny, arcL, curv)
%TRANSLATEPOLYPB1 Translates PolygonDiscretize data to solver data, and
%   produces boundary information.
%INPUTS
%   xq, vector NxM, The x coordinates of the quadrature points.
%   yq, vector NxM, The y coordinates of the quadrature points.
%   markQ, vector NxM, The mark to apply to the edge starting at each
%       point. 0 is no flow, 1 is high pressure, and 2 is low pressure.
%   nx, vector NxM, The x coordinates of the normals at each point.
%   ny, vector NxM, The y coordinates of the normals at each point.
%   arcL, vector NxM, The arc distance to the next quadrature point.
%   curv, vector NxM, The curvature at each quadrature point.
%OUTPUTS
%   jac, vector NxM, arc length scaled by angular distance.
%   dirVs, vector NxM, The Dirichlet conditions; ignore for indices that
%       are not Dirichlet.
%   neuVs, vector NxM, The Neuman conditions, ignore for indices that are
%       Dirichlet.
%   dirNodes, integer vector Q, List of indices that have Dirichlet
%       conditions.
   M = length(markQ);
   jac = arcL * (M / (2*pi));
   dirVs = zeros(M,1);
   neuVs = zeros(M,1);
   dirNodes = [];
   di = 1;
   M = length(markQ);
   for i = 1:M
       switch markQ(i)
           case 1
               dirNodes(di) = i;
               di = di + 1;
               dirVs(i) = 1;
           case 2
               dirNodes(di) = i;
               di = di + 1;
               dirVs(i) = 0;
               
           case 3
               dirNodes(di) = i;
               di = di + 1;
               dirVs(i) = 0.5;
           otherwise
               neuVs(i) = 0;
       end
   end
end