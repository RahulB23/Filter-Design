close all;
clearvars;
clc;

%specifications
m = 26;
q_m = floor(0.1*m);
r_m = m - 10*q_m;
BL = 25+1.7*q_m + 6.1*r_m;
BH = BL + 20;
tbw = 4*10^3;

% Band Edge specs
fs1 = BL*10^3-tbw;
fp1 = BL*10^3;
fp2 = BH*10^3;
fs2 = BH*10^3+tbw;
f_samp = 330e3;

wp1=fp1*2*pi/f_samp;
wp2=fp2*2*pi/f_samp;
ws2=fs2*2*pi/f_samp;
dlt = 0.15;

% Kaiser paramters
A = -20*log10(dlt);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end
d_omega_by_pi = tbw*2/f_samp;
d_omega=d_omega_by_pi*pi;
N_min = ceil((A-8) / (2*2.285*(d_omega)));          


addend = 8;  
n=N_min + addend;


w1=(ws1+wp1)/2;
w2=(ws2+wp2)/2;
bp_ideal = ideal_lp(w2,2*n+1) - ideal_lp(w1,2*n+1);


kaiser_win_shifted = (kaiser(2*n+1,beta)) 
num_arr = [-n:1:n];
kaiser_win = kaiser_win_shifted(num_arr+n+1);  
FIR_BandPass = bp_ideal .* kaiser_win;         

fvtool(FIR_BandPass,'Analysis','freq');  

[H,f] = freqz(FIR_BandPass,1, 1024*1024,f_samp);
plot(f,abs(H))
hold on;
title("Magnitude Response")
xlabel("Hz")
ylabel("|H(f)|")
xline(fp1,'--m');
xline(fp2,'--m');
yline(0.85,'--m');
yline(0.15,'--m');
grid

figure()
plot(FIR_BandPass,'-o')
title('Time-domain Sequence')
