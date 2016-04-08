function [ xq, yq, markQ, nx, ny, arcL, curv ] = PolygonDiscretize( xp, yp, markP, N )
%POLYGONDISCRETIZE This will create quadrature nodes on a polygon.
%INPUTS
%   xp, vector M, The x coordinates of the vertices of the polygon. Should
%       be in counterclockwise order. May be concave.
%   yp, vector M, The y coordinates of the vertices of the polygon.
%   markP, vector M, A marker to apply to the edge starting at each point.
%   N, integer, Number of quadrature points per edge.
%OUTPUTS
%   xq, vector NxM, The x coordinates of the quadrature points.
%   yq, vector NxM, The y coordinates of the quadrature points.
%   markQ, vector NxM, The mark to apply to the edge starting at each
%       point.
%   nx, vector NxM, The x coordinates of the normals at each point.
%   ny, vector NxM, The y coordinates of the normals at each point.
%   arcL, vector NxM, The arc distance to the next quadrature point.
%   curv, vector NxM, The curvature at each quadrature point.
    M = length(xp);
    [xq, yq, markQ, nx, ny] = GenQuadPts(xp, yp, markP, N);
    arcL = FigureArcLengths(xp, yp, xq, yq, N);
    curv = zeros(N*M,1);
end

function [arcL] = FigureArcLengths(xp, yp, xq, yq, N)
%FIGUREARCLENGTHS This will figure out the distances between quadrature
%   points.
%INPUTS
%   xp, vector M, The x coordinates of the vertices of the polygon. Should
%       be in counterclockwise order. May be concave.
%   yp, vector M, The y coordinates of the vertices of the polygon.
%   xq, vector NxM, The x coordinates of the quadrature points.
%   yq, vector NxM, The y coordinates of the quadrature points.
%   N, integer, Number of quadrature points per edge.
%OUTPUTS
%   arcL, vector NxM, The arc distance to the next quadrature point.
    M = length(xp);
    arcL = zeros(N*M,1);
    cI = 1;
    for i = 1:M
        for j = 1:(N-1)
            x1 = xq(cI);
            y1 = yq(cI);
            x2 = xq(cI+1);
            y2 = yq(cI+1);
            arcL(cI) = norm([x2-x1,y2-y1]);
            cI = cI + 1;
        end
        j = N;
        x1 = xq(cI);
        y1 = yq(cI);
        cxp = xp(wrapInd(i+1,M));
        cyp = yp(wrapInd(i+1,M));
        x2 = xq(wrapInd(cI+1,N*M));
        y2 = yq(wrapInd(cI+1,N*M));
        arcL(cI) = norm([cxp-x1,cyp-y1]) + norm([x2-cxp,y2-cyp]);
        cI = cI + 1;
    end
end

function [xq, yq, markQ, nx, ny] = GenQuadPts(xp, yp, markP, N)
%GENQUADPTS This will create quadrature nodes on a polygon.
%INPUTS
%   xp, vector M, The x coordinates of the vertices of the polygon. Should
%       be in counterclockwise order. May be concave.
%   yp, vector M, The y coordinates of the vertices of the polygon.
%   markP, vector M, A marker to apply to the edge starting at each point.
%   N, integer, Number of quadrature points per edge.
%OUTPUTS
%   xq, vector NxM, The x coordinates of the quadrature points.
%   yq, vector NxM, The y coordinates of the quadrature points.
%   markQ, vector NxM, The mark to apply to the edge starting at each
%       point.
%   nx, vector NxM, The x coordinates of the normals at each point.
%   ny, vector NxM, The y coordinates of the normals at each point.
    M = length(xp);
    xq = zeros(N*M,1);
    yq = zeros(N*M,1);
    markQ = zeros(N*M,1);
    nx = zeros(N*M,1);
    ny = zeros(N*M,1);
    for i = 1:M
        il = (i-1)*N+1;
        ih = i*N;
        x1 = xp(i);
        y1 = yp(i);
        x2 = xp(wrapInd(i+1,M));
        y2 = yp(wrapInd(i+1,M));
        cx = linspace(x1,x2,2*N+1);
        cy = linspace(y1,y2,2*N+1);
        xq(il:ih) = cx(2:2:(2*N+1));
        yq(il:ih) = cy(2:2:(2*N+1));
        markQ(il:ih) = markP(i);
        cnx = y2 - y1;
        cny = -(x2 - x1);
        ln = norm([cnx,cny]);
        cnx = cnx / ln;
        cny = cny / ln;
        nx(il:ih) = cnx;
        ny(il:ih) = cny;
    end
end

function wI = wrapInd(ind, Len)
%WRAPIND This will wrap an index into a vector.
%INPUTS
%   ind, integer, The index to wrap.
%   Len, integer, The length of the vector to wrap around.
%OUTPUTS
%   wI, integer, The wrapped index.
    wI = mod(ind-1,Len)+1;
end