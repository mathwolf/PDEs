function linear_spline_2_graph(mesh_points, data_points)

%	Math 550, Homework 2, Problem 1 (vii)
%   Solve the Galerkin FEM problem with Dirchlet-Neumann BCs

ALPHA = 0.0;
BETA = 2.0;                     % endpoints of domain

x = [ALPHA mesh_points BETA];   % vector containing mesh points for basis
points = length(x);       

% create load vector
% make matrices sparse so computations are more efficient
F = sparse(points-1,1);
for i = 1:points-2
    current_hat = @(z) hat(x(i), x(i+1), x(i+2), z);
    integrand = @(z) g(z)*current_hat(z);
    F(i) = integrate(integrand, x(i), x(i+2), 10);
end
current_hat = @(z) right_hat(points-1, points, z);
integrand = @(z) g(z)*current_hat(z);
F(points-1) = 9*cos(6) + integrate(integrand, x(points-1), x(points), 10);

% create stiffness matrix
% use symmetry when calculating
A = sparse(points-2, points-2);
for i = 1:points-2
   for j = i:points-2
      if abs(i-j) == 0
          hat_prime_i = @(z) hat_prime(x(i), x(i+1), x(i+2), z);
          integrand = @(z) sigma(z)*(hat_prime_i(z))^2;
          A(i,j) = integrate(integrand, x(i), x(i+2), 10);
      elseif abs(i-j) == 1
          hat_prime_i = @(z) hat_prime(x(i), x(i+1), x(i+2), z);
          hat_prime_j = @(z) hat_prime(x(j), x(j+1), x(j+2), z);
          integrand = @(z) sigma(z)*hat_prime_i(z)* hat_prime_j(z);
          A(i,j) = integrate(integrand, x(i), x(i+2), 10);
          A(j,i) = A(i,j);     % by symmetery
      end
   end
end
% Calculate matrix entries for boundary conditions
hat_prime_i = @(z) right_hat_prime(x(points-1), x(points), z);
integrand = @(z) sigma(z)*(hat_prime_i(z))^2;
A(points-1, points-1) = integrate(integrand, x(points-1), x(points), 10);

hat_prime_j = @(z) hat_prime(x(points-2), x(points-1), x(points), z);
integrand = @(z) sigma(z)*hat_prime_i(z)* hat_prime_j(z);
A(points-1, points-2) = integrate(integrand, x(points-1), x(points), 10);
A(points-2, points-1) = A(points-1, points-2); 

% Solve the system using Cholesky factorization
% After this step our approximate solution uh(z) below is well-defined
R = chol(A);
u = R'\F;
u = R\u;

% create approximate and exact data sets for plot
y_approx = zeros(length(data_points), 1);
for i = 1:length(y_approx)
   y_approx(i) = uh(data_points(i)); 
end

y_exact = zeros(length(data_points), 1);
for i = 1:length(y_approx)
   y_exact(i) = sin(pi*data_points(i));
end

% display the approximate solution
figure(1);
plot(data_points, y_approx);
title('FEM approximation to elliptic Dirichlet-Neumann problem');
xlabel('x');
ylabel('uh(x)');

% display error 
figure(2);
error = zeros(length(data_points), 1);
error = y_approx - y_exact;
plot(data_points, error);
title('Error in approximation to elliptic Dirichlet-Neumann problem');
xlabel('x');
ylabel('Error');



% Function containing the FEM approximation
    function y = uh(z)
       y = 0;
       for k = 1:points-2
           y = y + u(k) * hat(x(k), x(k+1), x(k+2), z);
       end
       y = y + u(points-1) * right_hat(x(points-1), x(points), z);
    end

end

