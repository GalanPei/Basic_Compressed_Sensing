function Exercise3_1
load('data3_1.mat');

plot(data3_1(1:40, 1), data3_1(1:40, 2)/100, 'o-')
xlabel('row number of matrix');
ylabel('recovery probabilities');
hold on;
plot(data3_1(41:80, 1), data3_1(41:80, 2)/100, 'o-')
legend('$k = 20$', '$k = 50$', 'interpreter', 'latex');
end