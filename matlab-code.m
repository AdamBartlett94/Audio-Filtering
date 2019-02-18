% Adv Digital Control - Assignment 1
% By Adam Bartlett
% ID: a1646071

%% Clear Workspace 
clear all
close all
clc


%% Initialisation
% Load initial recording
load('jimi2.mat')
initial_rec = x;
num_samples = size(x);

% Recording Parameters
Fs = 8000;          % Sampling frequency (Hz)
T = 1/Fs;           % Sampling period (s)
L = num_samples(1); % Length of signal
t = (0:L-1)*T;      % Time vector


%% Unfiltered Recording
% Listen to the unfiltered recording
display('Unfiltered Recording...')
sound(initial_rec, Fs);

% Plot: Unfiltered Recording: Time-Domain
figure
plot(t, initial_rec)
axis([0 5.1 -1 1])
title('Unfiltered Recording: Time-Domain')
xlabel('Time [s]')
ylabel('Sample [amp]')


%% Characterisation of Disturbance Signal
initial_rec_fft = fft(initial_rec);

% Using example from MATLAB's help information on fft function
P2 = abs(initial_rec_fft/L);
P1 = P2(1:(L/2+1));
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

% Plot: Unfiltered Recording: Single-Sided Amplitude Spectrum
figure
plot(f,P1) 
title('Unfiltered Recording: Single-Sided Amplitude Spectrum')
axis([0 4000 0 0.45])
xlabel('f (Hz)')
ylabel('|P1(f)|')

% Plot: Zoomed in version
figure
plot(f,P1) 
title('Unfiltered Recording: Single-Sided Amplitude Spectrum')
axis([110 130 0 0.45])
xlabel('f (Hz)')
ylabel('|P1(f)|')


%% Notch Filter: Transfer function determination
z = tf('z', T);
BW = 10;    % Bandwidth (Hz)
r = 1 - (pi/(2*4000/BW));

% Transfer Function:
H = ((z - cos(5.4*(pi/180)) - sin(5.4*(pi/180))*i)*(z - cos(5.4*(pi/180)) + sin(5.4*(pi/180))*i))/((z - r*cos(5.4*(pi/180)) - r*sin(5.4*(pi/180))*i)*(z - r*cos(5.4*(pi/180)) + r*sin(5.4*(pi/180))*i));

% Plot: Notch Filter - Bode Diagram
figure
bode(H)
title('Notch Filter - Bode Diagram')


%% Testing - Notch Filter
% Testing Parameters
L2 = 800;
t2 = (0:L2-1)*T;

% Generate 60Hz, 120Hz and 240Hz sin waves
wave_60 = sin(2*pi*60*t2);
wave_120 = sin(2*pi*120*t2);
wave_240 = sin(2*pi*240*t2);

%----------------------------Unfiltered-----------------------------------%
% Plot: Sine Wave Test - Unfiltered
figure
plot(t2, wave_60, t2, wave_120, t2, wave_240)
axis([0 inf -1.5 1.5])
legend('60Hz', '120Hz', '240Hz')
title('Sine Wave Test - Unfiltered')
xlabel('Time [s]')
ylabel('Sample [amp]')

%------------------------------Filtered-----------------------------------%
% Pass the sine waves through the notch filter
new_wave_60 = lsim(H, wave_60, t2);
new_wave_120 = lsim(H, wave_120, t2);
new_wave_240 = lsim(H, wave_240, t2);

% Plot: Sine Wave Test - Filtered
figure
plot(t2, new_wave_60, t2, new_wave_120, t2, new_wave_240)
axis([0 inf -1.5 1.5])
legend('60Hz', '120Hz', '240Hz')
title('Sine Wave Test - Filtered')
xlabel('Time [s]')
ylabel('Sample [amp]')


%% Filtering Recording
% Pass the initial recording through the notch filter
final_rec = lsim(H, initial_rec, t);

% Plot: Filtered Recording: Time-Domain
figure
plot(t, final_rec)
axis([0 5.1 -1 1])
title('Filtered Recording: Time-Domain')
xlabel('Time [s]')
ylabel('Sample [amp]')

% Listen to the filtered recording
pause(5)
display('Filtered Recording...')
sound(final_rec, 8000)
