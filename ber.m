clear all;
close all;
carrier_count = 10;
symbols_per_carrier = 3000;
train_size = 2000;
bits_per_symbol = 3;%16QAM
IFFT_bin_length = 512;
PrefixRatio = 1/4;
GI = PrefixRatio * IFFT_bin_length;
beta = 1/32;
GIP = beta*(IFFT_bin_length+GI);
SNR = 2;
tic;

%OFDM 信号产生
baseband_out_length = carrier_count * symbols_per_carrier * bits_per_symbol;%传输的总比特数
carriers = (1:carrier_count) + (floor(IFFT_bin_length/4) - floor(carrier_count/2));
conjugate_carriers = IFFT_bin_length - carriers + 2;
rand('twister',0);
baseband_out = round(rand(1,baseband_out_length));
complex_carrier_matrix = qam16(baseband_out);%16QAM
complex_carrier_matrix = reshape(complex_carrier_matrix',carrier_count,symbols_per_carrier)';
figure(1);
plot(complex_carrier_matrix,'*r');
axis([-4,4,-4,4]);
title('S8PSK星座图');
grid on;
IFFT_modulation = zeros(symbols_per_carrier,IFFT_bin_length);
IFFT_modulation(:,carriers) = complex_carrier_matrix;
IFFT_modulation(:,conjugate_carriers) = conj(complex_carrier_matrix);
signal_after_IFFT = ifft(IFFT_modulation,IFFT_bin_length,2);
time_wave_matrix = signal_after_IFFT;
figure(2);
plot(0:IFFT_bin_length-1,time_wave_matrix(2,:));
axis([0,512,-0.4,0.4]);
grid on;
ylabel('Amplitude');xlabel('time');
title('OFDM-S8PSK Time Signal, One Symbol Period');
XX = zeros(symbols_per_carrier,IFFT_bin_length+GI+GIP);
for k = 1:symbols_per_carrier
    for i = 1:IFFT_bin_length
        XX(k,i+GI) = signal_after_IFFT(k,i);
    end
    for i = 1:GI
        XX(k,i) = signal_after_IFFT(k,i+IFFT_bin_length-GI);
    end
    for j = 1:GIP
        XX(k,IFFT_bin_length+GI+j) = signal_after_IFFT(k,j);
    end
end
time_wave_matrix_cp = XX;
% figure(3);
% plot(0:length(time_wave_matrix_cp)-1,time_wave_matrix_cp(2,:));
% axis([0,600,-0.3,0.3]);
% grid on;
% ylabel('Amplitude');
% xlabel('Time');
% title('OFDM-S8PSK Time Signal with CP, One Symbol Period');
%加窗
windowed_time_wave_matrix_cp = zeros(1,IFFT_bin_length+GI+GIP);
for i = 1:symbols_per_carrier
    windowed_time_wave_matrix_cp(i,:) = real(time_wave_matrix_cp(i,:)).*rcoswindow(beta,IFFT_bin_length+GI)';
end
figure(4);
plot(0:IFFT_bin_length-1+GI+GIP,windowed_time_wave_matrix_cp(2,:));
axis([0,700,-0.2,0.2]);
grid on;
ylabel('Amplitude');
xlabel('Time');
title('OFDM-S8PSK Time Signal Apply a Window');
windowed_Tx_data = zeros(1,symbols_per_carrier*(IFFT_bin_length+GI)+GIP);
windowed_Tx_data(1:IFFT_bin_length+GI+GIP) = windowed_time_wave_matrix_cp(1,:);
for i = 1:symbols_per_carrier-1
    windowed_Tx_data((IFFT_bin_length+GI)*i+1:(IFFT_bin_length+GI)*(i+1)+GIP)...%
        = windowed_time_wave_matrix_cp(i+1,:);
end
Tx_data_withoutwindow = reshape(time_wave_matrix_cp',(symbols_per_carrier)...%
    *(IFFT_bin_length+GI+GIP),1)';
Tx_data = reshape(windowed_time_wave_matrix_cp',(symbols_per_carrier)*(IFFT_bin_length...%
    +GI+GIP),1)';
temp_time1 = (symbols_per_carrier)*(IFFT_bin_length+GI+GIP);
% figure(5);
% subplot(2,1,1);
% plot(0:temp_time1-1,Tx_data);
% grid on
% ylabel('Amplitude')
% xlabel('Time(samples)');
% title('OFDM-S8PSK Time Signal');
temp_time2 = symbols_per_carrier*(IFFT_bin_length+GI)+GIP;
subplot(2,1,2);
plot(0:temp_time2-1,windowed_Tx_data);
grid on;
ylabel('Amplitude')
xlabel('Time');
title('OFDM-S8PSK Time Signal');
symbols_per_average = ceil(symbols_per_carrier/5);
avg_temp_time = (IFFT_bin_length+GI+GIP)*symbols_per_average;
averages = floor(temp_time1/avg_temp_time);
average_fft(1:avg_temp_time) = 0;
for a = 0:(averages-1)
    subset_ofdm = Tx_data_withoutwindow(((a*avg_temp_time)+1):((a+1)*avg_temp_time));
    subset_ofdm_f = abs(fft(subset_ofdm));
    average_fft = average_fft+(subset_ofdm_f/averages);
end
average_fft_log = 20*log10(average_fft);
% figure(6);
% subplot(2,1,1);
% plot((0:(avg_temp_time-1))/avg_temp_time,average_fft_log);
% hold on;
% grid on;
% axis([0 0.5 -20 max(average_fft_log)]);
% ylabel('Magnitude(dB)');
% xlabel('Normalized Frequency(0.5 = fs/2)');
% title('OFDM-S8PSK Signal Spectrum');
symbols_per_average = ceil(symbols_per_carrier/5);
avg_temp_time = (IFFT_bin_length+GI+GIP)*symbols_per_average;
averages = floor(temp_time1/avg_temp_time);
average_fft(1:avg_temp_time) = 0;
for a = 0:(averages-1)
    subset_ofdm = Tx_data(((a*avg_temp_time)+1):((a+1)*avg_temp_time));
    subset_ofdm_f = abs(fft(subset_ofdm));
    average_fft = average_fft+(subset_ofdm_f/averages);
end
average_fft_log = 20*log10(average_fft);
subplot(2,1,2)
plot((0:(avg_temp_time-1))/avg_temp_time,average_fft_log);
hold on;
grid on;
axis([0 0.5 -20 max(average_fft_log)]);
ylabel('Magnitude(dB)');% figure(6);
% subplot(2,1,1);
% plot((0:(avg_temp_time-1))/avg_temp_time,average_fft_log);
% hold on;
% grid on;
% axis([0 0.5 -20 max(average_fft_log)]);
% ylabel('Magnitude(dB)');
% xlabel('Normalized Frequency(0.5 = fs/2)');
% title('OFDM-S8PSK Signal Spectrum');

xlabel('Windowed OFDM-S8PSK Signal Spectrum');
%NOISE
Tx_signal_power = var(windowed_Tx_data);
linear_SNR = 10^(SNR/10);
noise_sigma = Tx_signal_power/linear_SNR;
noise_scale_factor = sqrt(noise_sigma);
noise = randn(1,((symbols_per_carrier)*(IFFT_bin_length+GI))+GIP)*noise_scale_factor;
Rx_data = windowed_Tx_data+noise;
Rx_data_matrix = zeros(symbols_per_carrier,IFFT_bin_length+GI+GIP);
for i= 1:symbols_per_carrier
    Rx_data_matrix(i,:) = Rx_data(1,(i-1)*(IFFT_bin_length+GI)+1:i*(IFFT_bin_length...%
        +GI)+GIP);
end
Rx_data_complex_matrix = Rx_data_matrix(:,GI+1:GI+IFFT_bin_length);
Y1 = fft(Rx_data_complex_matrix,IFFT_bin_length,2);
Rx_carriers = Y1(:,carriers);
Rx_phase = angle(Rx_carriers);
Rx_mag = abs(Rx_carriers);
[M,N] = pol2cart(Rx_phase,Rx_mag);
Rx_complex_carrier_matrix = complex(M,N);
figure(7);
plot(Rx_complex_carrier_matrix,'*r');
axis([-4,4,-4,4]);
title('SNR = 30dB 接收数据星座图');
grid on;
Rx_serial_complex_symbols = reshape(Rx_complex_carrier_matrix',...%
    size(Rx_complex_carrier_matrix,1)*size(Rx_complex_carrier_matrix,2),1)';
amount = 0;

data.train_x = Rx_carriers(1:train_size,:);
data.test_x = Rx_carriers(train_size+1:symbols_per_carrier,:);
data.train_y = complex_carrier_matrix(1:train_size,:);
data.test_y = complex_carrier_matrix(train_size+1:symbols_per_carrier,:);
save('mydata_180902_2','data');

for i = 1:1:symbols_per_carrier
    for j = 1:1:carrier_count
        if(abs(Rx_carriers(i,j)-complex_carrier_matrix(i,j))>pi/8)
            amount = amount + 1;
        end
    end
end
% receive_out = zeros(1,baseband_out_length);
% m = 1;
% triples = [0,0,0;0,0,1;0,1,0;0,1,1;1,0,0;1,0,1;1,1,0;1,1,1];
% complex_array = [1.0, 0.7071+0.7071i, -0.7071+0.7071i, 1i, 0.7071-0.7071i, -1i, -1, -0.7071-0.7071i];
% results = zeros(symbols_per_carrier,carrier_count*3);
% inputs = reshape(baseband_out',carrier_count*3,symbols_per_carrier)';
% for i = 1:1:symbols_per_carrier
%     for j = 1:1:carrier_count
%         for k = 1:1:8
%             if(abs(complex_array(k)-Rx_carriers(i,j))<pi/8)
%                 results(i,(j-1)*3+1) = triples(k,1);
%                 results(i,(j-1)*3+2) = triples(k,2);
%                 results(i,(j-1)*3+3) = triples(k,3);
%                 break
%             end
%         end
%         m = m+1;
%     end
% end
toc;

