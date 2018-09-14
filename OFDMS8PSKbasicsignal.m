clear all;
close all;
P = 10;%number of carriers
M = 20;%number of symbols
T = 10e-6;
deltf = 1/T;
h = 1/4;
tx = 0:T/10:20*T;
signals_diff = [[0 -1 -2 -3 -4 3 2 1]; zeros(7,8)];
phi = zeros(M,P);
s = zeros(1,201);
%S8PSKÔ¤±àÂë
for i=2:8
    for j=1:8
        if j==1
            if(signals_diff(i-1,8) ~= -4)
                signals_diff(i,j)=signals_diff(i-1,8);
            end
            if(signals_diff(i-1,8) == -4)
                signals_diff(i,j)=4;
            end
        end
        if j~=1
            signals_diff(i,j)=signals_diff(i-1,j-1);
        end
    end
end

before = randsrc(1,M*P,[1:1:8]);
after = randsrc(1,M*P,[1:1:8]);
before = reshape(before,M,P);
after = reshape(after,M,P);
j = 0;
for t = 0:T/10:20*T
    j = j+1;
    for k = 1:1:M %number of symbols
        
        for p =1:1:P
            for q = 1:1:k
                phi(q,p) = phi(q,p) + pi*h*signals_diff(before(q,p),after(q,p));
            end
            if (t>= k*T && t<(k+1)*T)
                s(j) = s(j) + exp(sqrt(-1)*2*pi*p*t./T)*exp(sqrt(-1)*phi(k,p))/sqrt(T);
            end
            if (t<k*T || t>=(k+1)*T)
                    s(j) = s(j) + 0;
            end
        end
    end
end
plot(tx,s);
axis([0,2e-4,-2500,2500]);