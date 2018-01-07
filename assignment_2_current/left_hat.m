function [ y ] = left_hat( left, middle, x)
%   1d hat function on the left boundary
%   middle in this function is the boundary point

if x <= middle && x >= left
    y = (x - middle) / (left - middle);
else
    y = 0;
end

