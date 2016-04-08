function A = mixed_poisson_SLP_matrix(x,y,nx,ny,jac,cur,dirichlet_nodes)
%MIXED_POISSON_SLP_MATRIX This generates the matrix to solve.
%INPUTS
%	x, vector N, The x coordinates of the boundary.
%	y, vector N, The y coordinates of the boundary.
%	nx, vector N, The x coordinates of the normals at the boundary.
%	ny, vector N, The y coordinates of the normals at the boundary.
%	jac, vector N, ds/dt, arc-length scaled by the arc's angle.
%	cur, vector N, Curvature at each boundary node.
%	dirichlet_nodes, integer vector Q, List of node indices that have a Dirichlet condition.
%OUTPUTS
%	A, matrix NxN, The matrix to use to get the weight function.
	n = length(x);
	A = zeros(n,n);

	k = 6;

	switch k
		case 2
			gamma = [1.825748064736159e+00; -1.325748064736159e+00];
			
		case 6

			gamma = [4.967362978287758e+00, -1.620501504859126e+01, ...
						2.585153761832639e+01, -2.222599466791883e+01, ...
						9.930104998037539e+00, -1.817995878141594e+00]';
	end

	for i =  1:n
	   if (max(i == dirichlet_nodes))
		   for j = 1:n
			  if (i~=j)
				  dist = sqrt((x(i) - x(j))^2 + (y(i) - y(j))^2);
				  w = 0;
				  if (abs(i - j) <= k) %Kapur-Rokhlin rule
						w = gamma(abs(i-j));
				  end
				  
				  A(i,j) = (2*pi/(n))*(w+1)*log(dist)*jac(j)/(2*pi);              
			  end
		   end
	   else
		   for j = 1:n
			   
			   if (i == j)
				   A(i,j) = -(2*pi/n)*(cur(j)/2)*jac(j)/(2*pi) - 0.5;           
			   else
				   dist2 = (x(i) - x(j))^2 + (y(i) - y(j))^2;
				   rdotn = nx(i)*(x(i) - x(j)) + ny(i)*(y(i) - y(j));

				   A(i,j) = (2*pi/n)*rdotn*jac(j)/(dist2*2*pi);
			   end
		   end
	   end
	end


