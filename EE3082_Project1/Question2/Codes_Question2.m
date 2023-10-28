%% Create a message (2a) NOTE: I recommend that you run each style separately

Fs = 40000; % Sampling rate
td = 1/Fs; % Time Duration
T = 4; 
t = td:td:T; % approximation of cont time
% nBits = 16; % Number of bits to represent each sample
% nChannels = 1; % Mono channel
% ID = -1; % Default audio input device
% recObj1 = audiorecorder(Fs,nBits,nChannels,ID);
% disp('Start speaking first message for 4 seconds')
% recordblocking(recObj1,4);
% disp('End of recording the first message')
% myrecording = getaudiodata(recObj1);
% m1 = getaudiodata(recObj1);
% filename = 'message1.wav'; % Name the file 
% audiowrite(filename,m1,Fs);
[m1,Fs] = audioread("message1.wav"); % assign the message to m1
sound(m1,Fs)

% recObj2 = audiorecorder(Fs,nBits,nChannels,ID);
% disp('Start speaking second message for 4 seconds')
% recordblocking(recObj2,4);
% disp('End of recording the second message')
% myrecording = getaudiodata(recObj2);
% m2 = getaudiodata(recObj2);
% filename = 'message2.wav'; % Name the file 
% audiowrite(filename,m2,Fs);
[m2,Fs] = audioread("message2.wav"); % assign the message to m2
sound(m2,Fs)

% let's plot m1 and m2
figure(1)
subplot(211)
plot(t,m1)
title('Message1 Signal')
xlabel('Time-Seconds')
subplot(212)
plot(t,m2)
title('Message2 Signal')
xlabel('Time-Seconds')
%% Find Spectrum Analyzer 2(b)
Fc=8000; % carrier freq.
qam_mod = m1.*cos(2*pi*Fc*t')+m2.*sin(2*pi*Fc*t');
sa=dsp.SpectrumAnalyzer('SampleRate',Fs, ...
    'PlotAsTwoSidedSpectrum',true,'NumInputPorts',1, ...
    'ChannelNames',{'QAM Spectrum'});
sa(qam_mod);
release(sa);
%% Add Channel Noise (2c)
N = 1e-3; % it is valid for e
%N = 1e-6;
w = sqrt(N)*randn(size(qam_mod));
qam_mod = qam_mod + w;
% let's plot with m1
figure(2)
plot(t,m1)
hold on
plot(t,qam_mod)
title('Original Signal Modulated Signal with Noise')
xlim([2 2.01])
%% QAM demodulator and recover messages (2d)
m1_rec = 2*cos(2*pi*Fc*t').*qam_mod;
m2_rec = 2*sin(2*pi*Fc*t').*qam_mod;

% Design a LPF for Reconstruction
f_cutoff1=5000;
f_stop1=7000;
lpFilt1=designfilt('lowpassfir','PassbandFrequency',f_cutoff1,'StopbandFrequency',f_stop1,'samplerate',Fs);
fvtool(lpFilt1)

f_cutoff2 = 7000;
f_stop2 = 9000;
lpFilt2 = designfilt('lowpassfir','passbandfrequency',f_cutoff2,'stopbandfrequency',f_stop2,'samplerate',Fs);
fvtool(lpFilt2)

m1_demod_rec = filter(lpFilt1,m1_rec);
m2_demod_rec=filter(lpFilt2,m2_rec);

% Let's listen original messages and recovered messages
disp('Hit Enter to listen Original Message1')
pause
sound(m1,Fs)
disp('Hit Enter to listen Recovered Message1')
pause
sound(m1_demod_rec,Fs)
disp('Hit Enter to listen Original Message2')
pause
sound(m2,Fs)
disp('Hit Enter to listen Recovered Message2')
pause
sound(m2_demod_rec,Fs)

% let's plot 
figure(3)
plot(t,m1)
hold on
plot(t,m1_demod_rec)
legend('Original m1','Recovered m1','Location','best')
title('Original Message1 and Recovered Message1')
xlim([2 2.01])

figure(4)
plot(t,m2)
hold on
plot(t,m2_demod_rec)
legend('Original m2','Recovered m2','Location','best')
title('Original Message2 and Recovered Message2')
xlim([2 2.01])
