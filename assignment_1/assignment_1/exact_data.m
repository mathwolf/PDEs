function [u] = exact_data(A, B, N, T, M)
% Exact data used for analysis of Crank-Nicolson method
% Crank-Nicolson method for solution of the 1D parabolic problem
% Input parameters
%   A: lower endpoint of space domain
%   B: upper endpoint
%   N: number of mesh points in the spatial grid
%   T: upper end of time scale - here we assume time starts at 0
%   M: number of steps in time grid
% Output
%   Detailed explanation goes here

dt = T/M;                           % timestep
t = dt * [0:M];

h = (B - A)/N;                      % space step
x = A + h*[0:N]';

u = (sin(4*x) + cos(2*x)) * (cos(t) + sin(t));

end

