function [] = assignment_1_iii()
% Math 550, assignment 1
% Problem (iii)

k = 0.0;      % parameter for mesh
N = 10.0 * 2.0^k;
M = 2.0 * 10.0 * 2.0^k;     % twice as many for this problem

[x,t,u] = Crank_Nicolson('initial_data', 'source_term', ...
    'left_boundary_data', 'right_boundary_data', 0.1, 2, 3, N, 2, M);
surf(x,t,u');

end

