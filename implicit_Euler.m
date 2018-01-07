function [ x, t, U ] = implicit_Euler( vstring, fstring, a, A, B, N, T, M )
% Implicit Euler method for solution of the 1D parabolic problem
% Input parameters
%   vstring: initial condition
%   fstring: driving term
%   a: positive parameter in heat equation
%   A: lower endpoint of space domain
%   B: upper endpoint
%   N: number of mesh points in the spatial grid
%   T: upper end of time scale - here we assume time starts at 0
%   M: number of steps in time grid
% Output
%   Detailed explanation goes here

dt = 1/M;               % timestep
t = dt * [0:M];

h = (B - A)/N;          % space step
x = A + h*[0:N]';

U = zeros(N+1, M+1);    % zero BC is automatically filled in here
U(:,1) = feval(vstring(x));

iteration_matrix = speye(N-1) + (a*dt/h^2) * ...
    spdiags([-ones(N-1,1),2*ones(N-1,1),-ones(N-1,1)], [-1,0,1], N-1, N-1);
R = chol(iteration_matrix);

n = 2:N;
for m = 2:M+1
   b = U(n,m-1) + dt * feval(fstring, x(n), t(m));
   b = R'\b;
   U(n,m) = R\b;
end

end

