%% Initialisation
clc ; clear ; close all;
%% Listing all variable
load('data-proj.mat')
whos 
%% Plot of angular speed
fig=1;
figure(fig)
plot(t,omega)
grid on 
hold on
xlabel('Time [sec]')
ylabel('Angular speed [rad/sec]')
%% Sampling period
Te1=t(5)-t(4);
Te2= 0.05 ;%sampling period 
Fe2=1/Te2 ;%sampling frequency
Tf=t(end); %duration of the signal
N=Tf/Te2 ; %number of samples 
% (N=Duration/sampling period)

%frequency vector
f1=-Fe2*(N/2-1)/N:Fe2/N:0;
f2=Fe2/N:Fe2/N:(N/2)*Fe2/N;
f=[f2,f1];
O= zeros(N,1);
for m=1:N
  for k=1:N
    O(m)=O(m)+omega(k)*exp(-1i*2*pi*m*k/N);
  
  end
end
figure(fig+1);
stem(f,abs(O)/N)
grid on
xlim([-2,2])
xlabel('f [Hz]')
ylabel('DFT(\omega (t))')

% % Identify positive frequencies above a threshold
% threshold_amplitude = 0.2;
% positive_frequencies = f(abs(O) > threshold_amplitude);
% 
% % Deduce maximum frequency
% max_frequency = max(positive_frequencies);
% 
% fprintf('Maximum Frequency: %.2f Hz\n', max_frequency);

%create the new time vector 
t1=0:Te2:t(end)-Te2;
%% filter design
fc1=2;
H1=tf(1,[1/(2*pi*fc1)  1]);
Of=lsim(H1,O,t1);
%% plot of filtered signal
figure(1);
plot(t1,Of,'r')
grid on
legend(' \omega(t) unfiltered','\omega_{f}(t) filtered','Fontsize',14)

% dc_gain = dcgain(H1);
% fprintf('DC Gain of the Filter: %.4f\n', dc_gain);
% figure(fig+1);
% bode(H1);

