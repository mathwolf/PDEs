function [ I ] = integrate( f, a, b, n )
% Numerical estimation of integral in 1D using midpoint rule
% We subdivide the interval into n uniform pieces

% Input
% input_function    function to integrate, passed in as a string
% a                 left endpoint
% b                 right endpoint
% n                 number of intervals

% f = str2func(input_function);
I = 0;
length = (b-a) / n;
for i = 1:n
    midpoint = a + (i-1) * length + length/2;
    term = length * f(midpoint);
    I = I + term;
end

end

