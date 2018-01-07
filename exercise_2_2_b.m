function [] = exercise_2_2_b()
% Math 550, exercise sheet 2
% Problem 2, (iii)

[x,t,u] = explicit_Euler('initial_data', 'source_term', ...
    1.0, 1.0, 10, 1.0, 500);

figure
j = 1
for i = 1:2:11
   subplot(5,1,j);
   plot(x,u(:,i));
   j = j + 1;
end

end

