function [demodu_bit_symbol] = demodS8PSK(Rx_serial_complex_symbols)
N_sample = 4;
N = 1000;
h = 1/4;
m = 1;p = 4;
fc = 2;
Qn=generateQn(m,p);
Ts=1;%�źų�������
dt=Ts/N_sample;%�������
An=[4 2;1 3;2 4;3 1];%��ʾ��ǰ״̬��ת�ƾ���
out=zeros(4,N);out1=zeros(4,N);out2=zeros(4,N);%out1,out2Ϊת�Ƶ�ĳ��״̬ʱ��֧1�ͷ�֧2�ľ�����·����outΪת�Ƶ���Ӧ·��������С����������·��
sout=zeros(4,N+1);sout1=zeros(4,N+1);sout2=zeros(4,N+1);%sout1��sout2Ϊת�Ƶ�ĳ��״̬�ķ�֧����1�ͷ�֧����2���϶�Ӧ��·������������ֵ��soutΪ������ֵ�н�С��·������
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
      sout1(An(st,i+1),L+1)=sum+sout(st,L);%·������ֵ����ǰ��·������ֵ���Ϸ�֧����ֵ
      out1(An(st,i+1),1:L)=[out(st,1:L-1) 2*i-1];%��¼�µ�ǰ��·��
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
      sout2(An(st,i+1),L+1)=sum+sout(st,L);%·������ֵ����ǰ��·������ֵ���Ϸ�֧����ֵ
      out2(An(st,i+1),1:L)=[out(st,1:L-1) 2*i-1];%��¼�µ�ǰ��·��
      end
  end
    %�������ж�·����������������֧�����Ľϴ�ֵ����¼·��
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
a=findmax(A);%�ҳ�����ֵ
demodu_bit_symbol=out(a,:);
% end

