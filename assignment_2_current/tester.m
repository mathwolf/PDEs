function [] = tester()
%   Test FEM



ALPHA = 0.0;
BETA = 4.0;
x = [ALPHA 1. 2. 3. BETA];   % test mesh

points = length(x);             % number of mesh points

h = zeros(points-1, 1);         % length of each mesh element
for i = 1:length(h)
   h(i) = x(i+1) - x(i); 
end

% load vector
% make all matrices sparse
F = sparse(points-2,1);
for i = 1:points-2
    % Integrate using midpoint rule with one point
    current_hat = @(z) hat(x(i), x(i+1), x(i+2), z);
    integrand = @(z) f(z)*current_hat(z);
    F(i) = integrate(integrand, x(i), x(i+2), 5);
end

% stiffness matrix
% use symmetry and tridiagonal form when calculating
A = sparse(points-2, points-2);
for i = 1:points-2
   for j = i:points-2
      if abs(i-j) == 0
          % Integrate using midpoint rule with one point
          hat_prime_i = @(z) hat_prime(x(i), x(i+1), x(i+2), z);
          integrand = @(z) sigma(z)*(hat_prime_i(z))^2;
          A(i,j) = integrate(integrand, x(i), x(i+2), 50);
      elseif abs(i-j) == 1
          % Integrate using midpoint rule with one point
          hat_prime_i = @(z) hat_prime(x(i), x(i+1), x(i+2), z);
          hat_prime_j = @(z) hat_prime(x(j), x(j+1), x(j+2), z);
          integrand = @(z) sigma(z)*hat_prime_i(z)* hat_prime_j(z);
          A(i,j) = integrate(integrand, x(i), x(i+2), 50);
          A(j,i) = A(i,j);     % by symmetery
      end
   end
end

% Solve system using Cholesky
R = chol(A);
u = R'\F;
u = R\u;


disp(F);
disp(A);
disp(u);

% plot approximation
x_data = linspace(ALPHA, BETA, 100);
for i = 1:length(x_data)
   y_data = uh(x_data); 
end
plot(x_data,y_data);

% Nested function to compute approximation
% so that we do not need to pass the basis functions as arguments
    function y = uh(z)
        y = 0;
       for k = 1:points-2
           y = y + u(k) * hat(x(k), x(k+1), x(k+2), z);
       end
    end

end

