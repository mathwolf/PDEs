function [ y ] = hat( left, middle, right, x )
%	1d hat function
    if x <= middle
        if x >= left
            y = (x - left) / (middle - left);
        else
            y = 0;
        end
    else
        if x <= right
            y = (right - x) / (right - middle);
        else
            y = 0;
        end
    end
end

