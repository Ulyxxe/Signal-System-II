clear
close all
clc

%% Data visualization

% Question 1
load("data-proj.mat")
whos


%% Question 2
figure(1)
plot(t, omega)
grid on
hold on
xlabel('Time [sec]')
ylabel('Angular speed [rad/sec]')



%% Analog filtering

% Question 3
Te1 = t(2)-t(1);


%% Question 4
Te2= 0.05; 
Fe1=1/Te1;
Tf=t(end);
N=Tf/Te1;

f1=-Fe1*(N/2-1)/N:Fe1/N:0;
f2=Fe1/N:Fe1/N:(N/2)*Fe1/N;
f = [f2,f1];
w= zeros(N,1);
for m=1:N
  for k=1:N
    w(m)=w(m)+omega(k)*exp(-1i*2*pi*m*k/N);
  
  end
end

figure(2)
stem(f,abs(w)/N)
grid on
xlim([-2 2])
xlabel('f [Hz]')
ylabel('Amplitude Spectrum of Angular Speed Signal')


%% Question 6
% filter design
t1=0:Te1:t(end)-Te1;
fc=0.2;
wc=2*pi*fc;

H1=tf(1,[1/(2*pi*fc)  1]);
wf=lsim(H1,omega,t);

% plot of filtered signal
figure(1);
plot(t,wf,'r')
hold off
grid on
legend(' omega(t) unfiltered','omega_{f}(t) filtered','Fontsize',14)

data = [t',omega'];


%% Question 7
wf1=zeros(N,1);
for m = 1 : N
    for k = 1 : N
        wf1(m) = wf1(m) + wf(k) * exp(-1i*2*pi*m*k/N);
    end
end

figure(3)
subplot(2,1,1);
stem(f,abs(w)/N), hold on
stem(f,abs(wf1)/N, 'r'), hold off
grid on
xlim([-100 100])
legend(' DFT of omega(t)','DFT of omega_{f}(t)','Fontsize',14)
title('DFT(w(t) and DFT(w_{f}(t) within [-100,100] Hz')

subplot(2,1,2);
stem(f,abs(w)/N), hold on
stem(f,abs(wf1)/N, 'r'), hold off
grid on
xlim([-2 2])
legend(' DFT of omega(t)','DFT of omega_{f}(t)','Fontsize',14)
title('DFT(w(t) and DFT(w_{f}(t) within [-2,2] Hz')



%% Sampling

% Question 8 and 10
temp1 = 1:round(Te2/Te1):length(t);
we=wf(temp1);
Te = t(temp1);


%% Question 11
figure(4)
plot(t,wf), hold on
xlim([10 12])
xlabel('Time [sec]')
ylabel('Angular Velocity [rad/sec]')
grid on
stem(Te,we, 'r'), hold off
legend(' omega_{f}(t)','omega_{e}(t)','Fontsize',14)


%% Question 12
Fe2=1/Te2;
Tf2=Te(end);
N2=Tf2/Te2;

f3=-Fe2*(N2/2-1)/N2:Fe2/N2:0;
f4=Fe2/N2:Fe2/N2:(N2/2)*Fe2/N2;
f_2=[f4,f3];

we_dft = zeros(N2,1);
for m = 1 : N2
    for k = 1 : N2
        we_dft(m) = we_dft(m) + we(k) * exp(-1i*2*pi*m*k/N2);
    end
end

figure(5)
stem(f,abs(w)/N), hold on
stem(f,abs(wf1)/N, 'r')
stem(f_2,abs(we_dft)/N2, 'g'), hold off
xlabel('Frequency [Hz]')
ylabel('DFT')
legend({'DFT of omega(t)', 'DFT of omega_f(t)', 'DFT of omega_e(t)'})
xlim([-2 2])


%% Questions 13
% angular acceleration
wd_start=(we(2)-we(1))/Te2;
wd_end=(we(end)-we(end-1))/Te2;
wd_mid=zeros(N2-2,1);
for i=2:N2-1
    wd_mid(i)=(we(i+1)-we(i-1))/(2*Te2);
end

wd=[wd_start;wd_mid;wd_end];


% angular position
theta=zeros(N2,1);

for i=1:N2
    for k=1:i
        theta(i)=theta(i)+Te2*we(k);
    end
end


%% Question 14
t_ang=0:Te2:Te(end)-Te2;

figure(6)
plot(t_ang, theta)
xlabel('Time [sec]')
ylabel('Angular position [rad]')
grid on

figure(7)
plot(Te,wd)
xlabel('Time [sec]')
ylabel('Angular acceleration [rad/s^2]')
grid on


%% Question 15
wd_dft=zeros(N2,1);
for m=1:N2
    for k=1:N2
        wd_dft(m)=wd_dft(m)+wd(k)*exp(-1i*2*pi*m*k/N2);
    end
end

theta_dft = zeros(N2,1);
for m=1:N2
    for k=1:N2
        theta_dft(m)=theta_dft(m)+theta(k)*exp(-1i*2*pi*m*k/N2);
    end
end

figure(8)
stem(f_2,abs(theta_dft)/N2)
xlim([-2 2])
grid on
xlabel('Frequency [Hz]')
ylabel('DFT of theta(t)')

figure(9)
stem(f_2,abs(wd_dft)/N2)
xlim([-2 2])
grid on
xlabel('Frequency [Hz]')
ylabel("DFT of omega'(t)")


%% Question 16
[num,denum]=tfdata(H1,'v');
[num_digital, denum_digital] = bilinear(num,denum,Fe2,fc);

H2=tf(num_digital, denum_digital, Te2, 'Variable', 'z');


%% Question 17
wd_f=lsim(H2,wd,Te);

figure(10)
plot(Te,wd), hold on
plot(Te,wd_f, 'r'), hold off
grid on
xlabel('Time [sec]')
ylabel('Angular acceleration [rad/s^2]')
legend({"omega'(t) unfiltered", "omega'(t) filtered"})


%% Question 18
wd_f_dft=zeros(N2,1);
for i=1:N2
    for k=1:N2
        wd_f_dft(i)=wd_f_dft(i)+wd_f(k)*exp(-1i*2*pi*m*k/N2);
    end
end

figure(11)
stem(f_2,abs(wd_dft)/N2), hold on
stem(f_2,abs(wd_f_dft)/N2, 'r'), hold off
xlim([-2 2])
grid on
xlabel('Frequency [Hz]')
ylabel("DFT")
legend({"DFT of omega'(t) unfiltered", "DFT of omega'(t) filtered"})