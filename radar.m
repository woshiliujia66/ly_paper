clear all;
close all;
carrier_count = 20;%M
symbols_per_carrier = 20;%P
T = 10e-6;
deltf = 1/T;
for p = 0:1:symbols_per_carrier
    part1 = part1 + exp(-sqrt(-1)*2*pi*p*tau/T);
    for q = -(carrier_count-1):1:(carrier_count-1)
    for m = 0:1:M - abs(q)-1
        
    end
end