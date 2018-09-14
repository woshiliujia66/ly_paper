clear all;
close all;
x = [0,2,4,6,8,10,12,14,16];
y = [3225,2914,2383,1790,1093,519,146,27,1]/4000;
semilogy(x,y);
xlabel('SNR');
ylabel('SER');