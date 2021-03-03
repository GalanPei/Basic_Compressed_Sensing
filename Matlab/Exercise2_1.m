function Exercise2_1
load("data2_1.mat");
histogram(data2_1, 20, 'Normalization','pdf');
xlabel('$x$', 'interpreter', 'latex');
ylabel('$P(x)$', 'interpreter', 'latex');
hold on
y = -4:0.1:4;
f = exp(-y.^2./2)./(sqrt(2*pi));
plot(y,f,'black--','LineWidth',1.5)
legend('sample',...
'$y = \frac{1}{\sqrt{2\pi}}e^{-\frac{x^2}{2}}$', 'interpreter', 'latex',...
'fontsize', 16)
end