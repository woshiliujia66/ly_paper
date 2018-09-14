clear all;
close all;
P = 10;%number of carriers
M = 20;%number of symbols
T = 10e-6;%duration of pluse
deltf = 1/T;
tx = linspace(-T/2,T/2,10001);
s = zeros(1,10001);
d = randsrc(1,M);
j = 0;
for t = -T/2:T/10000:T/2
    j = j+1;
    for p = 1:1:P
        for m = 1:1:M
            s(j) = s(j) + d(m)*exp(sqrt(-1)*2*pi*P*deltf*t)/sqrt(P);
        end
    end
end
plot(tx,s);