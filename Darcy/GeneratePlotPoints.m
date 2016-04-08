function [ xt, yt, ins ] = GeneratePlotPoints( xq, yq, xsp, ysp )
%GENERATEPLOTPOINTS This will generate points inside a polygon.
%INPUTS
%   xq, vector M, X-coordinates of points of the polygon.
%   yq, vector M, Y-coordinates of points of the polygon.
%   xsp, integer, The number of grid points in the x direction.
%   ysp, integer, The number of grid points in the y direction.
%OUTPUTS
%   xt, matrix xsp x ysp, X-coordinates of the grid.
%   yt, matrix xsp x ysp, Y-coordinates of the grid.
%   ins, boolean matrix xsp x ysp, Whether each point is inside the grid.
    lowX = min(xq);
    higX = max(xq);
    lowY = min(yq);
    higY = max(yq);
    xvals = linspace(lowX, higX, xsp);
    yvals = linspace(lowY, higY, ysp);
    [xt, yt] = meshgrid(xvals, yvals);
    ins = zeros(xsp, ysp);
    for i = 1:xsp
        for j = 1:ysp
            cx = xt(i,j);
            cy = yt(i,j);
            [in,on] = inpolygon(cx,cy,xq,yq);
            if in
                if on
                    isIn = false;
                else
                    isIn = true;
                end
            else
                isIn = false;
            end
            ins(i,j) = isIn;
        end
    end
end

