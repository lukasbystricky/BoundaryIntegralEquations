function [ xt, yt ] = GenerateInteriorPts( xq, yq, N )
%GENERATEINTERIORPTS This will generate points inside a polygon.
%INPUTS
%   xq, vector M, X-coordinates of points of the polygon.
%   yq, vector M, Y-coordinates of points of the polygon.
%   N, integer, Number of points to generate.
%OUTPUTS
%   xt, vector N, X-coordinates inside the polygon.
%   yt, vector N, Y-coordinates inside the polygon.
    lowX = min(xq);
    higX = max(xq);
    lowY = min(yq);
    higY = max(yq);
    xt = zeros(N,1);
    yt = zeros(N,1);
    for i = 1:N
        cx = 0;
        cy = 0;
        isOut = true;
        while(isOut)
            r1 = rand();
            r2 = rand();
            cx = (r1*higX + (1-r1)*lowX);
            cy = (r2*higY + (1-r2)*lowY);
            [in,on] = inpolygon(cx,cy,xq,yq);
            if in
                if on
                    isOut = true;
                else
                    isOut = false;
                end
            else
                isOut = true;
            end
        end
        xt(i) = cx;
        yt(i) = cy;
    end
end

