function [  ] = hat_prime_tester(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

x = linspace(0., 10., 100);
y = zeros(100, 1);

current_hat_prime = @(z) hat_prime(1., 2., 3., z);

for i = 1:100
    y(i) = current_hat_prime(x(i));
end

plot(x,y);
end

