function y = g(x)
% Forcing term
if x < 0.5
    y = pi * cos(pi * x) + (1.5 - x) * pi^2 * sin(pi * x);
else
    y = -pi * cos(pi * x) + (0.5 + x) * pi^2 * sin(pi * x);
end

end

