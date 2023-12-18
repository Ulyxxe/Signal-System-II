%% initialization
clc
clear all
close all
%% time vector
t=0:0.001:10; %create the time vector
n=length(t); %calculate the length of I
tt=t'; %transpose of t
%% signal generation
s=sin(2*pi*t) + 0.5*cos(4*pi*t);
fig=1;
figure(fig)
plot(t,s)
grid on
xlabel('time(secs)')
ylabel('s(t)')
%%
Te=t(2)-t(1); %sampling period
Fe=1/Te; %sampling fequency
Tf=t(end);
N=Tf/Te;
%Now we define the frequency vector
f1=-Fe*(N/2-1)/N:Fe/N:0;
f2=Fe/N:Fe/N:(N/2)*Fe/N;
f=[f2,f1];
S=zeros(N,1);
for m=1:N
    for k=1:N
        S(m)=S(m)+s(k)*exp(-j*2*pi*m*k/N);
    end
end
%%
figure(fig+1)
stem(f,abs(S)/N)
grid
xlim([-3 3])
xlabel('freq(Hz)')
ylabel('DFT of s(t)')
%% add randim noise to s(t) to get sb(t)
vmin=-1;
vmax=+1;
sb=zeros(1,n);
for k=1:n
    sb(k)=s(k)+(vmax-vmin)*rand;
end
%% plot s(t) and sb(t)
fig=fig+1
figure(fig)
plot(t,s,'b')
hold on
plot(t,sb,'y')
hold off
grid
xlabel('freq(Hz)')
ylabel('s(t) and sb(t)')
%%  DFT
Sb=zeros(N,1);
for m=1:N
    for k=1:N
        Sb(m)=Sb(m)+sb(k)*exp(-j*2*pi*m*k/N);
    end
end
%%
fig=fig+1
figure(fig)
stem(f,abs(S)/N,'y')
hold on
stem(f,abs(Sb)/N,'b')
hold off
grid
xlim([-10 10])
xlabel('freq(Hz)')
ylabel('s(t) and sb(t)')
%% filyer design
fc1=3;
H1=tf([1],[1/(2*pi*fc1) 1]);
sf=lsim(H1,sb,t);

%% signal prepared for Simulink
A = [t',sb'];

% fig=fig+1
% figure(fig)
% % plot(t,sb,'y')
% % hold on
% % plot(out.tout,out.A_out,'k')
% % hold on
% % plot(t,s,'b')
% hold off
% grid
% xlabel('time(secs)')
% ylabel('sb(t) and sf(t)')
%%
fig=fig+1
figure(fig)
plot(t,sb,'y')
hold on
plot(t,sf,'k')
hold on
plot(t,s,'b')
hold off
grid
xlabel('time(secs)')
ylabel('sb(t) and sf(t)')
%%  DFT of the filtered signal
Sf=zeros(N,1);
for m=1:N
    for k=1:N
        Sf(m)=Sf(m)+Sf(k)*exp(-j*2*pi*m*k/N);
    end
end
%% plot of the spectrum of the filtered siganl
fig=fig+1;
figure(fig)
stem(f,abs(Sb)/N,'r')
hold on
stem(f,abs(Sf)/N,'b')
hold off
grid
xlim([-10 +10])
xlabel('freq(Hz)')
ylabel('Sb and sf')
legend('with noise','filtered','Location','northeast','FontSize',18)


%% Question 11 - Increasing the sampling interval to 0.05 sec

Te2 = 0.05;

%New time vector 
te = 0:Te2:Tf;

%increase in sampling period
m = Te2/Te;
%New number of samples
ne = length(te);

% Zero initialisation
se = zeros(1,ne);

for k = 0:ne-1
    se(k+1) = sf(m*k+1);
end

%% Question 12 - Plot of se(t) in the continuous time

%%
fig=fig+1
figure(fig)
plot(te,se,'k')
grid
xlabel('time(secs)')
ylabel('se(t)')

%% DFT of the subsampled signal
% to achieve this we need to define a new frequency vector
Fe2 = 1/ Te2;
N2 = Tf/Te2;
fe1=-Fe2*(N2/2-1)/N2:Fe2/N2:0;
fe2=Fe2/N2:Fe2/N2:(N2/2)*Fe2/N2;
fe=[fe2,fe1];

Se = zeros(N2,1);
for m=1:N2
    for k=1:N2
        Se(m)=Se(m)+se(k)*exp(-j*2*pi*m*k/N2);
    end
end

%% Plot the spectrum using stem

fig=fig+1;
figure(fig)
stem(fe,abs(Se)/N2,'k')
grid
xlim([-10 +10])
xlabel('freq(Hz)')
ylabel('Se')


%% Question 13 - Differentiation using the digital approach

ds(1) = (se(2) - se(1))/Te2;
ds(ne)= (se(end)- se(end))/Te2
for k =2: ne - 1
    ds(k) = (se(k)- se(k-1))/Te2;
end

%% Continuous plot of differentiation
fig=fig+1;
figure(fig)
plot(te,ds)
grid

ylabel('ds(t)')

%%
dS = zeros(N2,1);
for m=1:N2
    for k=1:N2
        dS(m)=dS(m)+ds(k)*exp(-j*2*pi*m*k/N2);
    end
end

%%
fig=fig+1
figure(fig)
stem(fe,abs(dS)/N2)

%% integration 
is = zeros(N2,1);
is(1)=Te2*se(1);
for k=2:ne
    is(k)= Te2*se(k) + is(k-1);
end

%% Plot
fig=fig+1;
figure(fig)
plot(te,is)
grid

ylabel('is(t)')

%% DFT 

iS = zeros(N2,1);
for m=1:N2
    for k=1:N2
        iS(m)=iS(m)+is(k)*exp(-j*2*pi*m*k/N2);
    end
end

%% plot
fig=fig+1;
figure(fig)
stem(fe,abs(iS/N2))
grid

ylabel('iS(t)')

%% IIR Filter Design

z = tf("z", Te2);
K=1;
a =(pi*fc1*K*Te2)/(1+pi*fc1*Te2);
b = (pi*fc1*Te2-1)/(1+pi*fc1*Te2);
Hz= tf([a a],[1,b], Te2);

dsf = lsim(Hz,ds,te);

%% Plot (time domain)
fig=fig+1;
figure(fig)
plot(te,dsf,"b")
grid
hold on 
plot(te,dsf,'r')
hold off


%% DFT
dSf = zeros(N2,1);
for m=1:N2
    for k=1:N2
        dSf(m)=dSf(m)+dsf(k)*exp(-j*2*pi*m*k/N2);
    end
end
%% plot
fig=fig+1;
figure(fig)
stem(fe,abs(dSf/N2),'b')
grid
hold on 
stem(fe,abs(dSf)/N2,'r')
hold off
ylabel('dSf(t)')








