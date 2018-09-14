clear all;
close all;
h = [1,1/8,1/4,1/2];
N = 30;
alpha = zeros(1,N);gama = zeros(4,N);
for n = 1:1:N%关于alpha 与gamma关系的申明
    if (alpha(n) == 4)
        gama(4,n) = 1;gama(3,n) = 1;gama(2,n) = 1;gama(1,n) = 1;
    end
    if (alpha(n) == -4)
        gama(4,n) = -1;gama(3,n) = -1;gama(2,n) = -1;gama(1,n) = -1;
    end
    if (alpha(n) == 3)
        gama(4,n) = 1;gama(3,n) = 1;gama(2,n) = -1;gama(1,n) = 1;
    end
    if (alpha(n) == -3)
        gama(4,n) = -1;gama(3,n) = -1;gama(2,n) = 1;gama(1,n) = -1;
    end
    if (alpha(n) == 2)
        gama(4,n) = 1;gama(3,n) = -1;gama(2,n) = 1;gama(1,n) = 1;
    end
    if (alpha(n) == -2)
        gama(4,n) = -1;gama(3,n) = 1;gama(2,n) = -1;gama(1,n) = -1;
    end
    if (alpha(n) == 1)
        gama(4,n) = 1;gama(3,n) = -1;gama(2,n) = 1;gama(1,n) = -1;
    end
    if (alpha(n) == -1)
        gama(4,n) = -1;gama(3,n) = 1;gama(2,n) = -1;gama(1,n) = 1;
    end
    if (alpha(n) == 0)
        gama(4,n) = 1;gama(3,n) = -1;gama(2,n) = -1;gama(1,n) = -1;
    end
end
b = zeros(4,N);
%关于b的定义
for l = 1:1:4
    for  m = 1:1:N
        b(l,m) = exp(sqrt(-1)*pi/8*gama(l,m));
    end
end

%关于c的定义
for l = 1:1:4
    C(l) = sin(2*h(l)*pi*t./2)/sin(h(l)*pi);
end
for l = 1:1:4
    CT(l) = sin(2*h(l)*pi*(T-t)./2)/sin(h(l)*pi);
end
for n = 1:1:N
    sym t
    z(n) = z(n) + int(r(t).* (CT(1)*C(2)*C(3)*C(4) + C(1)*CT(2)*C(3)*C(4)...%
        +CT(1)*CT(2)*C(3)*C(4))
end
for k =1:1:11
    
end
