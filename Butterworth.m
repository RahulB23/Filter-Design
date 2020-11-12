close all;
clearvars;
clc;


% specifications
m = 26;
q_m = floor(0.1*m);
r_m = m - 10*q_m;
BL = 25+1.7*q_m + 6.1*r_m;
BH = BL + 20;
tbw = 4*10^3;

% 2. Band Edge specifications
fs1 = BL*10^3-tbw;
fp1 = BL*10^3;
fp2 = BH*10^3;
fs2 = BH*10^3+tbw;
f_samp = 330e3;
ws1_by_pi=fs1*2/f_samp;
wp1_by_pi=fp1*2/f_samp;
wp2_by_pi=fp2*2/f_samp;
ws2_by_pi=fs2*2/f_samp;

% 3. Transformed Band Edge specs using Bilinear Transformation         
Ws1 = tan(fs1/f_samp*pi) ;
Wp1 = tan(fp1/f_samp*pi);
Wp2 = tan(fp2/f_samp*pi);
Ws2 = tan(fs2/f_samp*pi);
B=Wp2-Wp1;
W0=sqrt(Wp1*Wp2);

% 4. Frequency transformation 
WLs1=(Ws1^2-W0^2)/(B*Ws1);
WLp1=(Wp1^2-W0^2)/(B*Wp1);
WLp2=(Wp2^2-W0^2)/(B*Wp2);
WLs2=(Ws2^2-W0^2)/(B*Ws2);

% 5. Lowpass specs
Ws=min(abs(WLs1), WLs2)
Wp=WLp2

% 6. butterworth lowpass transfer function
dlt=0.15;
D1=1/(1-dlt)^2-1;
D2=1/dlt^2-1;
N = ceil(0.5*log(D2/D1)/log(Ws/Wp));
Wc_min = Wp/(D1^(1/(2*N)));
Wc_max = Ws/(D2^(1/(2*N)));
Wc=1.0715;                                   %%%%

% plot poles
sk=zeros(2*N,1);
for k=1:(2*N)
    sk(k)=Wc*cos(pi/2 + pi*(2*k-1)/(2*N)) + j*Wc*sin(pi/2 + pi*(2*k-1)/(2*N));
end
figure(1),plot(real(sk),imag(sk),'rX')  % poles
title("Poles of the analog lowpass transfer function")
xlabel("Re(s)")
ylabel("Im(s)")
axis equal
grid on
sk1 = abs(sk(1));
t = linspace(0,2*pi,100);
hold on
plot(sk1*cos(t),sk1*sin(t),'b-')  

% 7. bandpass transfer function
[num,den] = zp2tf([],[sk(1:N)],Wc^N);  
size(num)
size(den)
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        
analog_bpf(s) = analog_lpf((s*s + W0*W0)/(B*s));        
%coeffs of analog bpf
[ns, ds] = numden(analog_bpf(s)); 
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          
kd = ds(1);    
kn = ns(1);    
ds = ds/kd;
ns = ns/kn;
k=kn/kd;

% 8. Discrete bpf
discrete_bpf(z) = analog_bpf((z-1)/(z+1));              
[nz, dz] = numden(discrete_bpf(z));                                        
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              
kd = dz(1);                                             
k = dz(1);    
dz = dz/k;
nz = nz/k;
fvtool(nz,dz,'Analysis','freq');                        


[~,p,~]=tf2zp(nz,dz);
figure;
plot(real(p),imag(p),'rX');
title("Poles of the transfer function")
xlabel("Re(z)")
ylabel("Im(z)")
axis equal
grid on
t = linspace(0,2*pi,1000);
hold on
plot(cos(t),sin(t),'b-') 


%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, f_samp);
figure;
plot(f,abs(H),'LineWidth',1);
hold on;
title("Magnitude Response")
xlabel("Hz")
ylabel("|H(f)|")
xline(fs1,'--m');
xline(fp1,'--m');
xline(fp2,'--m');
xline(fs2,'--m');
yline(0.85,'--m');
yline(0.15,'--m');
grid





