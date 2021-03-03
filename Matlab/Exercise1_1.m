function Exercise1_1
load("data1_1.mat");
% Initialize the cofficient matrix and objective function
A = [1, 2; 2, 1; -1, 0];
b = [10; -1; 0];
fun = @(x, y)norm(A*[x; y] - b)^2;
% Plot the figure of objective function
x = -4 : 1/10 : 1;
y = 0 : 1/10 : 8;
z = zeros(round(length(y)), round(length(x)));
for i = 1 : round(length(x))
    for j = 1 : round(length(y))
        z(j, i) = fun(x(i), y(j));
    end
end
figure();
surf(x, y, z);
hold on
% Plot trajectory of points x_i
sol_x = data1_1(1, :);
sol_y = data1_1(2, :);
sol_z = zeros(1, length(sol_x));
for k = 1 : length(sol_x)
    sol_z(k) = fun(sol_x(k), sol_y(k));
end
plot3(sol_x, sol_y, sol_z, 'w');
xlabel('$x$', 'interpreter', 'latex');
ylabel('$y$', 'interpreter', 'latex');
zlabel('$\left \| Ax-b \right \|^2$', 'interpreter', 'latex');
% Plot the contour map
figure();
contour(x, y, z, [11 15 20 50 100 150 200], 'ShowText', 'on');
hold on
plot(sol_x, sol_y, 'r');
xlabel('$x$', 'interpreter', 'latex');
ylabel('$y$', 'interpreter', 'latex');

end