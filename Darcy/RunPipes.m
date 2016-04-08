Ns = [2,4,8,16,32,64];
Nhigh = 32;
xp = [0,1,1,3,3,22,23,23,22,20,19,19,20,20,21,22,22,21,16,16,15,15,12,12,11,11,7,7,9,9,7,7,10,10,6,6,4,4,3,3,1,1,0];
yp = [0,0,2,2,0,0,1,4,5,5,4,2,2,3,4,3,2,1,1,5,5,1,1,5,5,1,1,2,2,3,3,4,4,5,5,1,1,5,5,3,3,5,5];
markP = [0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1];
[testX,testY] = GenerateInteriorPts(xp,yp,2048);
[xOmega, yOmega, inside] = GeneratePlotPoints(xp,yp,100,100);

%figure the best solution
[xa, ya, markQ, nx, ny, arcL, curv] = PolygonDiscretize(xp,yp,markP,Nhigh);
[jaca, dirVs, neuVs, dirNodes] = TranslatePolyPb1(xa, ya, markQ, nx, ny, arcL, curv);
realSig = SolveProblem(xa, ya, nx, ny, jaca, curv, dirNodes, dirVs, neuVs);
%do weaker solutions
fprintf('N,Error\n');
solution = zeros(size(xOmega));

% for N = Ns
%    [xq, yq, markQ, nx, ny, arcL, curv] = PolygonDiscretize(xp,yp,markP,N);
%    [jac, dirVs, neuVs, dirNodes] = TranslatePolyPb1(xq, yq, markQ, nx, ny, arcL, curv);
%    curSig = SolveProblem(xq, yq, nx, ny, jac, curv, dirNodes, dirVs, neuVs);
%    
%    
%    %curErr = FigureError(xq, yq, jac, curSig, xa, ya, jaca, realSig, testX, testY);
% %    M = N * length(xp);
% %    fprintf('%d,%e\n',M,curErr);
% end

for i = 1:100
    for j = 1:100
        if (inside(i,j))
            solution(i,j) = evaluate_SLP(xa,ya, xOmega(i,j), yOmega(i,j), realSig, jaca);
        else
            solution(i,j) = nan;
        end
    end
end

surf(solution);