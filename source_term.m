function [f] = source_term(x,t)

f = sin(pi*x) * exp(-2*t);

end

