function [ ] = exercise_1_4( )
% Math 550 ex 1.4
% Make a table to compare derivative to difference operators
h = 0.2;
v = [0.0 0.2 0.4 0.6 0.8 1.0]';
f_prime = cos(v);
delta_plus = (sin(v+h)-sin(v)) / h;
delta_minus = (sin(v) - sin(v-h)) / h;
delta_bar = (1/2)*(delta_plus+delta_minus);

disp([v, f_prime, delta_plus, delta_minus, delta_bar]);

end

