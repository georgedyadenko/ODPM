ratio = 812*20;

load 'back_test\straddles\test_data.mat'
hst_r=zeros(size(hst,1)-ratio+1,1);
for i = ratio:1:size(hst, 1)
    hst_r(i-ratio+1) = hst(i, 4) / hst(i-ratio+1, 4);
end;
a=normmix(hst_r(ratio*22:ratio*25), 2, 2, 'b', 100)
pause
a=normmix(hst_r(ratio*24:ratio*26), 2, 2, 'b', 100)
pause
a=normmix(hst_r(ratio*25:size(hst_r,1)), 5, 2, 'b', 100)
