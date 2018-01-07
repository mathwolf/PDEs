function [ x, t, U ] = Crank_Nicolson( vstring, fstring, g1string, ...
                                        g2string, a, A, B, N, T, M )
% Crank-Nicolson method for solution of a 1D parabolic problem

% Input parameters
%   vstring: initial condition
%   fstring: driving term
%   g1string: boundary condition at lower endpoint
%   g2string: boundary condition at upper endpoint
%   a: positive parameter in heat equation
%   A: lower endpoint of space domain
%   B: upper endpoint
%   N: number of mesh points in the spatial grid
%   T: upper end of time scale - here we assume time starts at 0
%   M: number of steps in time grid

% Output
%   x: vector of space gridpoints
%   t: vector of time gridpoints
%   U: array of approximate solutions to the PDE at gridpoints

dt = T/M;                           % timestep
t = dt * [0:M];

h = (B - A)/N;                      % space step
x = A + h*[0:N]';

U = zeros(N+1, M+1);    
U(:,1) = feval(vstring, x);         % fill in initial condition

% Use sparse matrix setup to make code more efficient
Bh = (1/h^2) * spdiags([-ones(N-1,1),2*ones(N-1,1),-ones(N-1,1)], ...
                        [-1,0,1], N-1, N-1);
R = chol(speye(N-1) + (a*dt/2)*Bh);

n = 2:N;
for m = 2:M+1
   b = (speye(N-1) - (a*dt/2)*Bh) * U(n,m-1) + ...
       dt * feval(fstring, x(n), (1/2)* (t(m-1)+t(m)));
   b = R'\b;
   U(n,m) = R\b;
   U(1,m) = feval(g1string, t(m));
   U(N+1,m) = feval(g2string, t(m));
end

end

