function[Qn]=generateQn(m,p)%函数功能是产生在调制指数h=m/p的所有相位状态
 if rem(m,2)==0
     for i=1:p,
         Qn(i)=(i-1)*m*pi/p;
     end
 else rem(m,2)==1
     for i=1:2*p
         Qn(i)=(i-1)*m*pi/p;
     end
 end