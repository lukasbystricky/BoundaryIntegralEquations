Ns = [2,4,8,16,32,64,128,256];
Nhigh = 512;
xp = [0,1,1,0];
yp = [0,0,1,1];
markP = [1,0,2,0];
[testX,testY] = GenerateInteriorPts(xp,yp,2048);
%figure the best solution
[xa, ya, markQ, nx, ny, arcL, curv] = PolygonDiscretize(xp,yp,markP,Nhigh);
[jaca, dirVs, neuVs, dirNodes] = TranslatePolyPb1(xa, ya, markQ, nx, ny, arcL, curv);
realSig = SolveProblem(xa, ya, nx, ny, jaca, curv, dirNodes, dirVs, neuVs);
%do weaker solutions
fprintf('N,Error\n');

xOmega = linspace(0.25,0.75);
yOmega = linspace(0.25,0.75);

solution = zeros(length(xOmega),length(yOmega));
err_trap = zeros(length(Ns),1);
f_handle = @(x,y) x.^2 - y.^2;
g_handle = @(x,y) 2*x.^2 - 2*y.^2;


for k = 1:length(Ns)
   [xq, yq, markQ, nx, ny, arcL, curv] = PolygonDiscretize(xp,yp,markP,Ns(k));
   [jac,dirVs, neuVs, dirNodes] = TranslatePolyPb1(xq, yq, markQ, nx, ny, arcL, curv);
   
%    dirNodes = 1:4*Ns(k);
%    [dirVs, neuVs] = create_bc(xq,yq, f_handle, g_handle);
   
   curSig = SolveProblem(xq, yq, nx, ny, jac, curv, dirNodes, dirVs, neuVs);
   
   for i = 1:length(xOmega)
       for j = 1:length(yOmega)
            solution(i,j) = evaluate_SLP(xq,yq,xOmega(i),yOmega(j),curSig,jac);
            err_trap(k) = err_trap(k) + (solution(i,j) - (1 - yOmega(j)))^2;%evaluate_SLP(xa,ya,xOmega(i),yOmega(j),realSig,jaca))^2;         
       end
   end
   
   err_trap(k) = sqrt(err_trap(k));
%    curErr = FigureError(xq, yq, jac, curSig, xa, ya, jaca, realSig, testX, testY);
%    M = N * length(xp);
%    fprintf('%d,%e\n',M,curErr);
end


figure();

surf(xOmega,yOmega,solution);
view(0,90)
shading interp

rates = zeros(length(Ns) -1,1);

for i = 1:length(rates)
   rates(i) = log(err_trap(i)/err_trap(i+1))/log(2); 
end