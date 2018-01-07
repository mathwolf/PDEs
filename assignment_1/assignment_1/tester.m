function [] = tester()
%   Test FEM



alpha = 0.0;
beta = 1.0;
x = [alpha 0.2 0.5 0.8 beta];
h = zeros(length(x)-1, 1);
for i = 1:length(h)
   h(i) = x(i+1) - x(i); 
end
n = length(h);

% load vector
F = zeros(n,1);
for i = 2:n+1
    % Integrate using midpoint rule with one point
    midpoint1 = (x(i) + x(i-1))/2;
    midpoint2 = (x(i+1) + x(i))/2;
    F(i-1) =  h(i) * f(midpoint1) * ...
            hat(x(i-1), x(i), x(i+1), midpoint1) + ...
            h(i+1) * f(midpoint2) * ...
            hat(x(i-1), x(i), x(i+1), midpoint2);
end

% stiffness matrix
A = zeros(n,n);
for i = 2:n+1
   for j = 2:n+1
      if abs(i-j) == 0
          % Integrate using midpoint rule with one point
          A(i-1,j-1) = (x(i) - x(i-1)) * sigma( (x(i) + x(i-1))/2 ) * ...
                hat_prime(x(i-1), x(i), x(i+1), (x(i) + x(i-1))/2 )^2 + ...
                (x(i+1) - x(i)) * sigma( (x(i+1) + x(i))/2 ) * ...
                hat_prime(x(i-1), x(i), x(i+1), (x(i+1) + x(i))/2 )^2;
      elseif abs(i-j) == 1
          % Integrate using midpoint rule with one point
          A(i-1,j-1) = (x(i) - x(i-1)) * sigma( (x(i) + x(i-1))/2 ) * ...
                hat_prime(x(i-1), x(i), x(i+1), (x(i) + x(i-1))/2 ) * ...
                hat_prime(x(j-1), x(j), x(j+1), (x(i) + x(i-1))/2 ) + ...
                (x(i+1) - x(i)) * sigma( (x(i+1) + x(i))/2 ) * ...
                hat_prime(x(i-1), x(i), x(i+1), (x(i+1) + x(i))/2 ) * ...
                hat_prime(x(j-1), x(j), x(j+1), (x(i+1) + x(i))/2 );

      end
   end
end

disp(F);
disp(A);
%disp(y);
%plot(x,y);

end

