function Exercise1_2

A_1 = [1 2; 2 1; -1 0];
A_2 = [1 2; 2 1; 1.8 -2];
A_3 = [1 2; 2 1; -2 -2];
A_4 = [1 2; 2 1; 2 -4];
A_5 = [1 2; 2 1; 3 -3];
cond_num = [cond(A_2'*A_2, 2) cond(A_5'*A_5, 2) cond(A_4'*A_4, 2) ...
            cond(A_1'*A_1, 2) cond(A_3'*A_3, 2)]
iter_num = [2 12 22 38 74];
plot(cond_num, iter_num)