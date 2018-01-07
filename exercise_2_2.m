function [] = exercise_2_2()
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

[x,t,u] = explicit_Euler('initial_data', 'source_term', ...
    1.0, 1.0, 10, 1.0, 500);
surf(x,t,u');
end

