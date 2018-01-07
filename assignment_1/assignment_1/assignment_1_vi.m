function [] = assignment_1_vi()
% Math 550, assignment 1
% Problem (vi)

% Store table data
table_data = zeros(6,4);

for j = 1:6
    
    k = j - 1.0;      % parameter for mesh
    N = 10.0 * 2.0^k;
    M = 2.0 * 10.0 * 2.0^k;     % twice as many for this dimension

    [x,t,U] = Crank_Nicolson('initial_data', 'source_term', ...
        'left_boundary_data', 'right_boundary_data', 0.1, 2, 3, N, 2, M);
    
    % Calculate  maximum error at final time
    exact = exact_data(2,3,N,2,M);
    max_err = max( abs( U(:,end)- exact(:,end) ) );
    % Store data for display in table
    table_data(j,1) = 0.1 / 2.0^(k);       % h 
    table_data(j,2) = 0.1 / 2.0^(k);       % delta t
    table_data(j,3) = max_err;              % e infinity, delta t, h
    table_data(j,4) = 2 * (0.1 / 2.0^(k))^2;    % order of convergence
 
    % Store calculated values for 3D plots if k=0 or k=4
    if k == 0
        x0 = x;
        t0 = t;
        U0 = U;
        exact0 = exact;
    elseif k == 4
        x4 = x;
        t4 = t;
        U4 = U;
        exact4 = exact;
    end
end

% display table of error data
disp('    Delta t     h       max err     EOC');
disp(table_data);

% Make 2D plot of error relationship
figure(1);
plot(table_data(:,2), table_data(:,3))
title('Error Relationship for Crank-Nicolson Method')
xlabel('Spatial grid size (h)')
ylabel('Maximum error at final time')

% Make the sets of 3D plots
figure(2);
colormap(gray);
surf(x0, t0, U0');
title('Crank-Nicolson solution for grid size h = 1/10')
xlabel('Spatial coordinate (x)')
ylabel('Temporal coordinate (t)')
zlabel('Calculated u(x,t)')

figure(3);
colormap(gray);
surf(x0, t0, exact0');
title('Exact solution for grid size h = 1/10')
xlabel('Spatial coordinate (x)')
ylabel('Temporal coordinate (t)')
zlabel('Exact u(x,t)')

figure(4);
colormap(gray);
surf(x0, t0, exact0' - U0');
title('Error for grid size h = 1/10')
xlabel('Spatial coordinate (x)')
ylabel('Temporal coordinate (t)')
zlabel('Exact value minus calculated')

figure(5);
colormap(gray);
surf(x4, t4, U4');
title('Crank-Nicolson solution for grid size h = 1/160')
xlabel('Spatial coordinate (x)')
ylabel('Temporal coordinate (t)')
zlabel('Calculated u(x,t)')

figure(6);
colormap(gray);
surf(x4, t4, exact4');
title('Exact solution for grid size h = 1/160')
xlabel('Spatial coordinate (x)')
ylabel('Temporal coordinate (t)')
zlabel('Exact u(x,t)')

figure(7);
colormap(gray);
surf(x4, t4, exact4' - U4');
title('Error for grid size h = 1/160')
xlabel('Spatial coordinate (x)')
ylabel('Temporal coordinate (t)')
zlabel('Exact solution minus calculated')

end

