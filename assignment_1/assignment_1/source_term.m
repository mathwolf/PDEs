function [f] = source_term(x, t)

ut = (sin(4.0*x) + cos(2.0*x)) * (-sin(t) + cos(t));
uxx = (-16.0*sin(4.0*x) - 4.0*cos(2.0*x)) * (cos(t) + sin(t));
f = ut - 0.1 * uxx;

end

