Ns = [2,4,8,16,32,64,128,256];
Nhigh = 512;
xp = [0,0,1,1];
yp = [0,1,1,0];
markP = [2,0,1,0];
[testX,testY] = GenerateInteriorPts(xp,yp,2048);
%figure the best solution
[xa, ya, markQ, nx, ny, arcL, curv] = PolygonDiscretize(xp,yp,markP,Nhigh);
[jaca, dirVs, neuVs, dirNodes] = TranslatePolyPb1(xa, ya, markQ, nx, ny, arcL, curv);
realSig = SolveProblem(xa, ya, nx, ny, jaca, curv, dirNodes, dirVs, neuVs);
%do weaker solutions
fprintf('N,Error\n');

xOmega = linspace(-20,20,100);
yOmega = linspace(-20,20,100);

pressure = zeros(length(xOmega),length(yOmega));
vel_u = zeros(length(xOmega), length(yOmega));
vel_v = zeros(length(xOmega), length(yOmega));

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
           if ((xOmega(i)  > 1 || xOmega(i) <  0) || (yOmega(j)  > 1 || yOmega(j) <  0))
            [pressure(i,j), vel_u(i,j), vel_v(i,j)] = evaluate_SLP(xq,yq,xOmega(i),yOmega(j),curSig,jac);
            
            err_trap(k) = err_trap(k) + (pressure(i,j) - (1 - yOmega(j)))^2;%evaluate_SLP(xa,ya,xOmega(i),yOmega(j),realSig,jaca))^2;  
           else
               pressure(i,j) = nan;
               vel_u(i,j) = nan;
               vel_v(i,j) = nan;
           end
       end
   end
   
   err_trap(k) = sqrt(err_trap(k));
%    curErr = FigureError(xq, yq, jac, curSig, xa, ya, jaca, realSig, testX, testY);
%    M = N * length(xp);
%    fprintf('%d,%e\n',M,curErr);
end


figure();

surf(xOmega,yOmega,vel_u);
view(0,90)
shading interp

rates = zeros(length(Ns) -1,1);

for i = 1:length(rates)
   rates(i) = log(err_trap(i)/err_trap(i+1))/log(2); 
end