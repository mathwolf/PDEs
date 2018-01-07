function [ y ] = right_hat( left, middle, x)
    %   1d hat function on the right boundary
    %   middle in this function is the boundary point
    if x >= left && x <= middle
        y = (x - left) / (middle - left);
    else
        y = 0; 
end

