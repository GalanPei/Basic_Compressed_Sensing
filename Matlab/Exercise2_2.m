function Exercise2_2
load('data2_2.mat');

for i = 1 : length(data2_2(:, 1))
    semilogy(1:length(data2_2(i, :)), data2_2(i, :))
    hold on
end
xlabel('Iteration step')
ylabel('Error')

end