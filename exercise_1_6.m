function [  ] = exercise_1_6(  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
h = [0.32 0.16 0.08 0.04 0.02 0.01]';

eh = zeros(6,1);
e2h = zeros(6,1);
e4h = zeros(6,1);
for j = 1:6
    eh(j,1) = (sin(1.0 + h(j,1)) - sin(1.0)) / h(j,1) - cos(1.0);
    e2h(j,1) = (sin(1.0 + 2*h(j,1)) - sin(1.0)) / (2*h(j,1)) - cos(1.0);
    e4h(j,1) = (sin(1.0 + 4*h(j,1)) - sin(1.0)) / (4*h(j,1)) - cos(1.0);
end
disp([h, eh, (eh/e2h), (e4h-e2h)/(e2h-eh)]);
% disp([h', eh', (eh/e2h)', ((e4h-e2h)/(e2h-eh))']);
end

