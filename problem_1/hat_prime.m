function [ y ] = hat_prime( left, middle, right, x )
%	Derivative of 1d hat function
    if x <= middle
        if x >= left
            y = 1 / (middle - left);
        else
            y = 0;
        end
    else
        if x <= right
            y = -1 / (right - middle);
        else
            y = 0;
        end
    end
end

