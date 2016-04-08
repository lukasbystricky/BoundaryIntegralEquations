% N = 100;
% theta = (0:N-1)'*2*pi/N;
% 
% x = cos(theta);
% y = sin(theta);
% 
% nx = cos(theta);
% ny = sin(theta);
% jac = ones(length(nx),1);
% cur = ones(length(ny),1);
% 
% dirichlet_nodes = 1:5;
% neumann_nodes = 6:N;
% 
% A = mixed_poisson_SLP_matrix(x,y,nx,ny,jac,cur,dirichlet_nodes);
% 
% I = zeros(size(A));
% I(neumann_nodes, neumann_nodes) = -0.5*eye(length(neumann_nodes));
% A = A + I;
% 
% f = @(x,y) 5*x + y;
% g = @(x,y) 5*x + y;
% 
% rhs(dirichlet_nodes) = f(x(dirichlet_nodes),y(dirichlet_nodes));
% rhs(neumann_nodes) = g(x(neumann_nodes), y(neumann_nodes));
% 
% sigma = A\rhs';
% 
% x1 = linspace(-0.5,0.5);
% y1 = linspace(-0.5,0.5);
% 
% p = zeros(length(x1),length(y1));
% 
% for i = 1:length(x1)
%     for j = 1:length(y1)
%         p(i,j) = evaluate_SLP(x, y, x1(i), y1(j), sigma, jac);
%     end
% end
% 
% surf(x1,y1,p);

%% calculate convergence

exact_solution = @(x,y) x.^2 - y.^2;
f = @(x,y) x.^2 - y.^2;
g = @(x,y) 2*x.^2 - 2*y.^2;

xOmega = linspace(-1.4,1.4);
yOmega = linspace(-1.4,1.4);

N = [128, 256, 512, 1024];
err_trap = zeros(length(N),1);
p_approx = zeros(length(xOmega),length(yOmega));
conditions = zeros(length(N),1);

for i = 1:length(N)
   
    rhs = zeros(N(i),1);
    theta = (0:N(i)-1)'*2*pi/N(i);
    
    r = 1 + 0.4*cos(5*theta);
    rp = -2*sin(5*theta);
    rpp = -10*cos(5*theta);
    
    
%     xBoundary = r.*cos(theta);
%     yBoundary = r.*sin(theta);
%     
%     nx = (rp.*sin(theta) + r.*cos(theta))./sqrt(rp.^2+r.^2);
%     ny = (r.*sin(theta) - rp.*cos(theta))./sqrt(rp.^2+r.^2);
%     jac = sqrt(rp.^2 + r.^2);
%     cur_n = r.^2 + 2*rp.^2 - r.*rpp;
%     cur_d = (r.^2 + rp.^2).^(3/2);
%     cur = cur_n./cur_d;

    xBoundary = cos(theta);
    yBoundary = sin(theta);
% 
    nx = cos(theta);
    ny = sin(theta);
    jac = ones(length(nx),1);
    cur = ones(length(ny),1);

   dirichlet_nodes = 1:N(i)/2; 
   neumann_nodes = N(i)/2 + 1:N(i);
%     neumann_nodes = [];
%    dirichlet_nodes = 1:N(i);
   
    
    A = mixed_poisson_SLP_matrix(xBoundary,yBoundary,nx,ny,jac,cur,dirichlet_nodes);
    
    rhs(dirichlet_nodes) = f(xBoundary(dirichlet_nodes),yBoundary(dirichlet_nodes));
    rhs(neumann_nodes) = g(xBoundary(neumann_nodes), yBoundary(neumann_nodes));
    
    sigma = (A)\rhs;    
    conditions(i) = cond(A);
    
    points = 0;
    
    for j = 1:length(xOmega)
        for k = 1:length(yOmega)
            %test if point is inside domain
            
            thetaTmp = atan2(yOmega(k),xOmega(j));
            xTmp1 = (0.5 + 0.4*cos(5*thetaTmp))*cos(thetaTmp);
            yTmp1 = (0.5 + 0.4*cos(5*thetaTmp))*sin(thetaTmp);
            
            xTmp2 = (0.5 + 0.4*cos(5*thetaTmp+pi))*cos(thetaTmp+pi);
            yTmp2 = (0.5 + 0.4*cos(5*thetaTmp+pi))*sin(thetaTmp+pi);
            
            if (xOmega(j) > min(xTmp1,xTmp2) && xOmega(j) < max(xTmp1,xTmp2) && ...
                    yOmega(k) > min(yTmp1,yTmp2) && yOmega(k) < max(yTmp1,yTmp2))
                
            %if (sqrt(xOmega(j)^2 + yOmega(k)^2)<0.5)
                points = points + 1;
                p_approx(j,k) = evaluate_SLP(xBoundary, yBoundary, xOmega(j), yOmega(k), sigma, jac);
                err_trap(i) = err_trap(i) +  ((exact_solution(xOmega(j),yOmega(k))- p_approx(j,k)))^2;
            end
        end
    end
    
    p_approx(p_approx == 0) = nan;    
    err_trap(i) = sqrt(err_trap(i)/points);
    
end

figure();

surf(xOmega,yOmega,p_approx);
view(0,90)
shading interp

rates = zeros(length(N) -1,1);

for i = 1:length(rates)
   rates(i) = log(err_trap(i)/err_trap(i+1))/log(2); 
end