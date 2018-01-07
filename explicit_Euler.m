function [ x, t, U ] = explicit_Euler( vstring, fstring, a, L, N, T, M)
% Math 550, exercise sheet 1
% Explicit Euler method for PDEs
% Input
%       vstring: initial data
%       fstring: driving term
%       a: parameter for parabolic PDE
%       L: total length
%       N: number of steps in spatial grid
%       T: total time
%       M: number of steps in time grid
% Output
%       U: matrix of calculated data point
%   Detailed explanation goes here

dt = T/M;                           % time step
t = dt * [0:M];                     % we start at t = 0

h = L/N;                            % spatial step
x = h * [0:N]';                     % also starting at 0

U = zeros(N+1, M+1);                % this fills in our 0 BCs
U(:,1) = feval(vstring, x);         % initial conditions

alpha = a * dt / (h^2);

n = 2:N;
for m = 1:M
   U(n,m+1) = U(n,m) + dt * feval(fstring, x(n), t(m)) + ...
       alpha * (U(n+1,m) - 2*U(n,m) + U(n-1,m));
end

end

