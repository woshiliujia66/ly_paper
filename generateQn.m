function[Qn]=generateQn(m,p)%���������ǲ����ڵ���ָ��h=m/p��������λ״̬
 if rem(m,2)==0
     for i=1:p,
         Qn(i)=(i-1)*m*pi/p;
     end
 else rem(m,2)==1
     for i=1:2*p
         Qn(i)=(i-1)*m*pi/p;
     end
 end