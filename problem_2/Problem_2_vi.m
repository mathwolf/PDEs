function Problem_2_vi()
% Set up a uniform mesh covering (0,2) with spacing h = 0.0025
mesh = linspace(0.0025, 1.9975, 799);
% Evaluate approximation at 1000 equally spaced points
evaluation = linspace(0, 1, 5000);

linear_spline_2_graph(mesh, evaluation);

end

