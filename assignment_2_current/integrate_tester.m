function [  ] = integrate_tester(  )

g = @(x) hat(0., 1., 2., x);
y = integrate(g, 0.0, 1.0, 500);
disp(y);
%disp(g(1));
%disp(g(1.5));
end

