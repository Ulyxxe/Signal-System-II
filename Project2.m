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
xlabel('Time [sec]')
ylabel('Angular speed [rad/sec]')
%% Sampling period
Te1=t(5)-t(4);
Te2= 0.05 ;%sampling period 
Fe2=1/Te2 ;%sampling frequency
Tf=t(end); %duration of the signal
N=Tf/Te2 ;    %number of samples (Duration/sampling period)
%frequency vector
f1=-Fe2*(N/2-1)/N:Fe2/N:0;
f2=Fe2/N:Fe2/N:(N/2)*Fe2/N;
f=[f2,f1];
O= zeros(1,N);
for m=1:N
  for k=1:N
    O(m)=O(m)+omega(k)*exp(-j*2*pi*m*k/N);
  end
end
fig = 2;
stem(f,abs(O)/N)
xlim([-2,2])
xlabel('f [Hz]')
ylabel('\omega (t)')