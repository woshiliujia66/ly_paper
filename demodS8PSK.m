function [demodu_bit_symbol] = demodS8PSK(Rx_serial_complex_symbols)
N_sample = 4;
N = 1000;
h = 1/4;
m = 1;p = 4;
fc = 2;
Qn=generateQn(m,p);
Ts=1;%信号持续周期
dt=Ts/N_sample;%抽样间隔
An=[4 2;1 3;2 4;3 1];%表示当前状态的转移矩阵
out=zeros(4,N);out1=zeros(4,N);out2=zeros(4,N);%out1,out2为转移到某个状态时分支1和分支2的经历的路径，out为转移到对应路径度量较小的所经历的路径
sout=zeros(4,N+1);sout1=zeros(4,N+1);sout2=zeros(4,N+1);%sout1，sout2为转移到某个状态的分支度量1和分支度量2加上对应的路径度量的两个值，sout为这两个值中较小的路径度量
 for L=1:N
     for st=1:2
     for i=0:1,
        x=(2*i-1)*ones(1,N_sample);
        t=L-1:dt:L-dt;
        s_phase(1)=Qn(st)/(pi*h);
           for m=1:N_sample-1,
           s_phase(m+1)=s_phase(m)+x(m)*dt;
           end
        s_cpfsk=cos(2*pi*fc*t+pi*h*s_phase);
        sum=0;
        for n=N_sample*(L-1)+1:N_sample*L,
            s=n-N_sample*(L-1);
            y = real(Rx_serial_complex_symbols(n));
            sum=sum+real(Rx_serial_complex_symbols(n))*s_cpfsk(s)*dt;
        end;
      sout1(An(st,i+1),L+1)=sum+sout(st,L);%路径度量值等于前的路径度量值加上分支度量值
      out1(An(st,i+1),1:L)=[out(st,1:L-1) 2*i-1];%记录下当前的路径
      end
  end
  for st=3:4
     for i=0:1,
        x=(2*i-1)*ones(1,N_sample);
        t=L-1:dt:L-dt;
        s_phase(1)=Qn(st)/(pi*h);
           for m=1:N_sample-1,
           s_phase(m+1)=s_phase(m)+x(m)*dt;
           end
        s_cpfsk=cos(2*pi*fc*t+pi*h*s_phase);
        sum=0;
        for n=N_sample*(L-1)+1:N_sample*L,
            s=n-N_sample*(L-1);
            sum=sum+real(Rx_serial_complex_symbols(n)).*s_cpfsk(s)*dt;
        end;
      sout2(An(st,i+1),L+1)=sum+sout(st,L);%路径度量值等于前的路径度量值加上分支度量值
      out2(An(st,i+1),1:L)=[out(st,1:L-1) 2*i-1];%记录下当前的路径
      end
  end
    %下面是判断路径度量加上两个分支度量的较大值并记录路径
            for i=1:4
                if sout1(i,L+1)>=sout2(i,L+1),
                    sout(i,L+1)=sout1(i,L+1);
                    out(i,1:L)=out1(i,1:L);
                else    
                    sout(i,L+1)=sout2(i,L+1);
                    out(i,1:L)=out2(i,1:L);
                end
            end
end
A=[sout(1,L+1) sout(2,L+1) sout(3,L+1 ) sout(4,L+1)];
a=findmax(A);%找出最大的值
demodu_bit_symbol=out(a,:);
% end

