function Problem_1_vii()
% Set up a uniform mesh with h = 0.01
mesh = linspace(0.01, 0.99, 99);
% Evaluate approximation at 5000 equally spaced points
evaluation = linspace(0, 1, 5000);

linear_spline_1_graph(mesh, evaluation);

end

