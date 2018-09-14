function [complex_qam_data] = qam16(bitdata)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
h = 1/4;
signals_diff = [[0 -1 -2 -3 -4 3 2 1]; zeros(7,8)];
phi = zeros(1,length(bitdata)/3);
%S8PSK预编码
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

if (bitdata == 0)
    bitdata = -1;
end
X1 = reshape(bitdata,3,length(bitdata)/3)';
alpha = zeros(1,length(bitdata)/3);
for i = 1:1:length(bitdata)/3
    X1(i,1) = X1(i,1);
    X1(i,2) = xor(X1(i,1),X1(i,2));
    X1(i,3) = xor(X1(i,2),X1(i,3));
    alpha(i) = 2^2*X1(i,1)+2^1*X1(i,2)+X1(i,1)+1;
end

phi(1) = 0;
for i =2:1:length(bitdata)/3
    for q = 1:1:i
        phi(i) = phi(i) + pi*h*signals_diff(alpha(i-1),alpha(i));
    end
end

for i = 1:1:length(bitdata)/3
    complex_qam_data(i) = exp(sqrt(-1)*phi(i));
end
% X1 = reshape(bitdata,4,length(bitdata)/4)';
% d = 1;
% for i = 1:length(bitdata)/4;
%     for j = 1:4
%         X1(i,j) = X1(i,j)*(2^(4-j));
%     end
%     source(i,1) = 1+sum(X1(i,:));
% end
% mapping = [-3*d 3*d; -d 3*d; d 3*d;3*d 3*d;-3*d d;-d d;d d;3*d d;-3*d -d;-d -d;d -d;3*d -d;
%     -3*d -3*d;-d -3*d;d -3*d;3*d -3*d];
% for i = 1:length(bitdata)/4
%     qam_data(i,:) = mapping(source(i),:);
% end
% complex_qam_data = complex(qam_data(:,1),qam_data(:,2));
