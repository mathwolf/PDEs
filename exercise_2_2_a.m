function [] = exercise_2_2_a()
% Math 550, exercise sheet 2
% Problem 2, (ii)

[x,t,u] = explicit_Euler('initial_data', 'source_term', ...
    1.0, 1.0, 10, 1.0, 500);
surf(x,t,u');
end

