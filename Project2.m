%% Initialisation
clc ; clear ; close all;
%% Listing all variable
load('data-proj.mat')
whos 
%% Plot of angular speed
fig=1
figure(fig)
plot(omega,t)
grid on 
xlabel("Time [sec]")
ylabel("Angular speed [rad/sec]")