%% Create a message (3a)

Fs = 40000; % Sampling rate
td = 1/Fs; % Time Duration
T = 4; 
t = td:td:T; % approximation of cont time
% nBits = 16; % Number of bits to represent each sample
% nChannels = 1; % Mono channel
% ID = -1; % Default audio input device
% recObj = audiorecorder(Fs,nBits,nChannels,ID);
% disp('Start speaking first message for 4 seconds')
% recordblocking(recObj,4);
% disp('End of recording the first message')
% myrecording = getaudiodata(recObj);
% m1 = getaudiodata(recObj);
% filename = 'message.wav'; % Name the file 
% audiowrite(filename,m1,Fs);
[m1,Fs] = audioread("message.wav"); % assign the message to m1
sound(m1,Fs);

% let's plot m1
figure(1)
plot(t,m1)
title('Message Signal')
xlabel('Time-Seconds')
%% Show the USB and DSB spectrum analyzer (3b)
Fc=8000; % carrier freq.
m1_hil = imag(hilbert(m1));

m1_usb_mod = m1.*cos(2*pi*Fc*t') - m1_hil.*sin(2*pi*Fc*t'); % USB modulation 
m1_dsb_mod = m1.*2.*cos(2*pi*Fc*t'); % DSB modulation

sa = dsp.SpectrumAnalyzer('SampleRate',Fs, ...
    'PlotAsTwoSidedSpectrum',true,'NumInputPorts',3, ...
    'ChannelNames',{'Message','USB signal','DSB signal'});
sa(m1,m1_usb_mod,m1_dsb_mod);
release(sa);
%% Add Channel Noise (3c)
 N = 1e-3; % it is valid for e
%N = 1e-6;
w = sqrt(N)*randn(size(m1));
usb_w = m1_usb_mod + w;
dsb_w = m1_dsb_mod + w;

% let's plot usb_w and dsb_w with m1
figure(2)
subplot(211)
plot(t,m1)
hold on
plot(t,usb_w)
xlim([1.5 1.503])
legend('Message Signal','USB Mod. Signal')
xlabel('Time-Seconds')
hold off

subplot(212)
plot(t,m1)
hold on
plot(t,dsb_w)
xlim([1.5 1.503])
legend('Message Signal','DSB Mod. Signal')
xlabel('Time-Seconds')
%% USB/LSB demod. and recover the message (3d)
usb_dem = 2*cos(2*pi*Fc*t').*usb_w;
lsb_dem = 2*cos(2*pi*Fc*t').*dsb_w;

% Design a LPF for Reconstruction
f_cutoff=5000;
f_stop=7000;
lpFilt=designfilt('lowpassfir','PassbandFrequency',f_cutoff,'StopbandFrequency',f_stop,'samplerate',Fs);
fvtool(lpFilt)

usb_m1_rec1 = filter(lpFilt,usb_dem); % Recover the m1 with usb_dem
lsb_m1_rec2 = filter(lpFilt,lsb_dem); % Recover the m1 with lsb_dem

disp('Hit Enter to listen Original Message')
pause
sound(m1,Fs)
disp('Hit Enter to listen Recovered Message with usb_dem')
pause
sound(usb_m1_rec1,Fs)
disp('Hit Enter to listen Recovered Message with lsb_dem')
pause
sound(lsb_m1_rec2,Fs)

% let's plot 
figure(3)
plot(t,m1)
hold on
plot(t,usb_m1_rec1)
legend('Original m1','Recovered m1 usb','Location','best')
title('Original Message and Recovered Message with USB')
xlim([1.5 1.503])

figure(4)
plot(t,m1)
hold on
plot(t,lsb_m1_rec2)
legend('Original m1','Recovered m1 lsb','Location','best')
title('Original Message and Recovered Message with LSB')
xlim([1.5 1.503])
