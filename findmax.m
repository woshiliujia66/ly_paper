function[n]=findmax(shulie)
%������һ��һά���飬������ҳ����һά�����еĵڼ����������ģ�������ڼ��������������������������
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
