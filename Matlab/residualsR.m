function residualsR
load("residuals.mat");
load("kvalues.mat")
for k = 0:3
subplot(2,2,k + 1);
for i = 30*k + 1:30*(k + 1)
    semilogy(1:1002, kvalues(i,:));
    hold on;
end
title("k = "+(-10*k + 40));
xlabel("iterations");
ylabel("residuals"); 
end
end

