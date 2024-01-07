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

%% QUESTION 7:
clear, close all, clc
% Parameters
Fe = 100; % Sampling frequency
Te = 1 / Fe; % Sampling period
tf = 8; % Total time
N = tf / Te; % Number of samples
t = 0 : Te : tf - Te; % Time vector
% Signal ω(t)
omega_t = cos(2 * pi * t) + cos(pi * t / 2);
% DFT of ω(t)
omega_t_dft = fft(omega_t);
omega_t_dft_shifted = fftshift(omega_t_dft);
freqs = (-N/2 : N/2 - 1) * (Fe / N); % Frequency vector
% Design a low-pass filter (Butterworth filter)
cutoff_freq = 5; % Cutoff frequency in Hz
order = 2; % Order of the filter
[b, a] = butter(order, cutoff_freq / (0.5 * Fe), 'low');
% Apply the filter to ω(t) to obtain ω_f(t)
omega_f_t = filtfilt(b, a, omega_t);
% DFT of ω_f(t)
omega_f_t_dft = fft(omega_f_t);
omega_f_t_dft_shifted = fftshift(omega_f_t_dft);
% Plot the amplitude spectrum of both ω(t) and ω_f(t) [-100 to +100 Hz]
figure;
stem(freqs, abs(omega_t_dft_shifted), 'b', 'filled');
hold on;
stem(freqs, abs(omega_f_t_dft_shifted), 'r', 'filled');
xlim([-100, 100]); % Limit frequency axis
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Amplitude Spectrum of ω(t) and ω_f(t) [-100 to +100 Hz]');
legend('ω(t)', 'ω_f(t)');
hold off;
% Plot the amplitude spectrum of both ω(t) and ω_f(t) [-2 to +2 Hz]
figure;
stem(freqs, abs(omega_t_dft_shifted), 'b', 'filled');
hold on;
stem(freqs, abs(omega_f_t_dft_shifted), 'r', 'filled');
xlim([-2, 2]); % Limit frequency axis
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Amplitude Spectrum of ω(t) and ω_f(t) [-2 to +2 Hz]');
legend('ω(t)', 'ω_f(t)');
hold off;

%% QUESTION 8/9:
% Sampling period Te2
Te2 = 0.05;
% Create a time vector for ωe(t)
t_e = 0:Te2:tf - Te2;
% Resample ωf(t) at the specified time instants
omega_e_t = interp1(t, omega_f_t, t_e);
size_omega_e_t = length(omega_e_t);
% Display the size of ωe(t)
disp(size_omega_e_t);
% Plot the resampled ωe(t)
figure;
plot(t_e, omega_e_t, 'b');
xlabel('Time (s)');
ylabel('Angular Velocity');
title('Resampled ωe(t)');

%% QUESTION 10:
% Sampling period Te2
Te2 = 0.05;
% Create a time vector for ωe(t)
t_e = 0:Te2:tf - Te2;
% Resample ωf(t) at the specified time instants
omega_e_t = interp1(t, omega_f_t, t_e);
% Create the new time vector te
te = 0:Te2:(length(omega_e_t) - 1) * Te2;
% Plot the resampled ωe(t) with the new time vector te
figure;
plot(te, omega_e_t, 'b');
xlabel('Time (s)');
ylabel('Angular Velocity');
title('Resampled ωe(t) with te');

%% QUESTION 11:
% Sampling period Te2
Te2 = 0.05;
% Total time
tf = 8; % Replace this with your actual total time
% Assuming t and omega_f_t are defined earlier in the code or passed as inputs
% Create a time vector for ωe(t)
t_e = 0:Te2:(tf - Te2);
% Resample ωf(t) at the specified time instants
omega_e_t = interp1(t, omega_f_t, t_e);
% Create a figure
figure;
% Plot both signals on the same graph with different colors
plot(t, omega_f_t, 'b', 'LineWidth', 1.5);  % Filtered angular velocity ωf(t) in blue
hold on;
plot(t_e, omega_e_t, 'r', 'LineWidth', 1.5);  % Sampled filtered angular velocity ωe(t) in red
xlabel('Time (s)');
ylabel('Angular Velocity');
title('Filtered Angular Velocity ωf(t) and Sampled ωe(t)');
legend('ωf(t)', 'ωe(t)');
grid on;
% Adjust xlim if necessary to fit within the range of your data
% xlim([10, 12]);
% Hold off to stop superimposing future plots
hold off;

%% QUESTION 12:
% Assuming omega_e_t, omega_t, and omega_f_t are defined
% Calculate the amplitude spectrum and frequency vector for omega_e_t
N_omega_e = length(omega_e_t);
f_omega_e = linspace(-0.5, 0.5, N_omega_e) * (1/Te2);
amp_spectrum_omega_e = abs(fftshift(fft(omega_e_t))) / N_omega_e;
% Calculate the amplitude spectrum and frequency vector for omega_t
N_omega = length(omega_t);
f_omega = linspace(-0.5, 0.5, N_omega) * (1/Te2);
amp_spectrum_omega = abs(fftshift(fft(omega_t))) / N_omega;
% Calculate the amplitude spectrum and frequency vector for omega_f_t
N_omega_f = length(omega_f_t);
f_omega_f = linspace(-0.5, 0.5, N_omega_f) * (1/Te2);
amp_spectrum_omega_f = abs(fftshift(fft(omega_f_t))) / N_omega_f;
% Plot the amplitude spectra
figure;
plot(f_omega_e, amp_spectrum_omega_e, 'r', 'LineWidth', 1.5); hold on;
plot(f_omega, amp_spectrum_omega, 'b', 'LineWidth', 1.5);
plot(f_omega_f, amp_spectrum_omega_f, 'g', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Amplitude Spectra of omega_e_t, omega_t, and omega_f_t');
legend('omega_e_t', 'omega_t', 'omega_f_t');
grid on;
hold off;

