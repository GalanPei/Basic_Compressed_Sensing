function Exercise2_3
load('data2_3.mat')
load('data2_4.mat')

figure()
plot(1 : length(data2_3), data2_3, 'o')
xlabel('$k$', 'interpreter', 'latex');
ylabel('Computional overhead');

end