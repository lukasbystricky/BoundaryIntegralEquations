function [ curErr ] = FigureError( xq, yq, jac, curSig, xa, ya, jaca, realSig, testX, testY )
%FIGUREERROR This will calculate the error.
%INPUTS
%   xq, vector N, The x coordinates of the approximate boundary.
%   yq, vector N, The y coordinates of the approximate boundary.
%   jac, vector N, The arclength scaled by the angular distance for the
%       approximate boundary.
%   curSig, vector N, The weight function on the approximate boundary.
%   xa, vector M, The x coordinates of the high density boundary.
%   ya, vector M, The y coordinates of the high density boundary.
%   jaca, vector M, The arclength scaled by the angular distance for the
%       high density boundary.
%   realSig, vector M, The weight function on the high density boundary.
%   testX, vector P, The x coordinates of the interior positions to test
%       at.
%   testY, vector P, The y coordinates of the interior positions to test
%       at.
%OUTPUTS
%   curErr, real, The error.
    
    %TODO
    
end

