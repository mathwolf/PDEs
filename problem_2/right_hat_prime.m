function y = right_hat_prime(left, middle, x)
% Derivative of 1d hat function on the right boundary
    if x >= left && x <= middle
        y = 1 / (middle - left);
    else
        y = 0; 
end

