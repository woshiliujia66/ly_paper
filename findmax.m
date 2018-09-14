function[n]=findmax(shulie)
%输入是一个一维数组，输出是找出这个一维数组中的第几个数是最大的，是输出第几个数，而不是输出这个最大数。
sym p;
Q=length(shulie);
max=shulie(1);
for i=1:Q
    if shulie(i)>=max,
        max=shulie(i);
        p=i;
    end
end
n=p;
