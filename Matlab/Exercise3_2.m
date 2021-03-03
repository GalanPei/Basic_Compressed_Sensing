function Exercise3_2
load('data3_2.mat');

imagesc([0.04 0.99], [0.04 0.99], 1/100*data3_2')
set(gca,'YDir','normal')
xlabel('Sparsity k/n');
ylabel('Number of equations m/n');
colorbar
end
