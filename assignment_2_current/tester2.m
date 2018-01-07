function [  ] = tester2(  )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

n=100;
x = linspace(0.0, 1.0, n);
y = zeros(n,1);
for i = 1:n
   y(i) = hat_prime(0., 0.5, 1., x(i)) 
end

plot(x,y);
end

