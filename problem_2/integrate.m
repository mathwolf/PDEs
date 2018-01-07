function [ I ] = integrate( f, a, b, n )
% Numerical integration in 1D using midpoint rule
% Subdivide the interval into n uniform pieces

% Input
% input_function    function to integrate
% a                 left endpoint
% b                 right endpoint
% n                 number of intervals

I = 0;
length = (b-a) / n;
for i = 1:n
    midpoint = a + (i-1) * length + length/2;
    term = length * f(midpoint);
    I = I + term;
end

end

