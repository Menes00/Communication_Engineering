%% Create a message (1a)  

Fs = 40000; % Sampling rate
td = 1/Fs; % Time Duration
T = 4; 
t = td:td:T; % approximation of cont time
nBits = 16; % Number of bits to represent each sample
nChannels = 1; % Mono channel
ID = -1; % Default audio input device
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

% plot message signal
figure(1)
plot(t,m1)
title('Message Signal')
xlabel('Time-Seconds')

%% Take a Downsample (1b)
Fs = 40000; % Sampling rate
td = 1/Fs; % Time Duration
T = 4; % 
t = td:td:T; % approximation of cont time
n = length(m1); % message length

fs = 8000; % we sample the(approximately) cont time signal at a rate fs. every 5 samples for 8kHz
ts = 1/fs;  % Sampling period
N = ts/td; % N should be an integer N = 5
s_out = downsample(m1,N);
s_out = upsample(s_out,N); 
% let's plot
figure(2)
plot(t,m1)
hold on
stem(t,s_out,'r')
grid on
xlabel('Time')
legend('Original Signal','Sampled Signal','Location','SouthWest')
title('Original and Sampled Signal')
xlim([2 2.001]) % very short time
%% find spectrum analyzer (1c)
sa = dsp.SpectrumAnalyzer('SampleRate',Fs,'PlotAsTwoSidedSpectrum',true, ...
    'NumInputPorts',2,'ChannelNames', ...
    {'Message Signal','Sampled Signal'});
sa(m1,s_out);
release(sa);
%% find low pass filter 1(d)
% take the Fourier Transform of m1 an s_out
f = (-(n-1)/2:(n-1)/2)*(Fs/n); % Generate the discrete frequency vector x axis 
fre_m1 = fftshift(fft(m1,n)); % Computes the FT of original (approximately) cont. time signal
fre_s_out = fftshift(fft(s_out,n)); % Computes the FT of sampled disc. time signal
% let's plot
figure(3)
stem(f,(fre_m1)/n,'b-s')
hold on 
stem(f,(fre_s_out)/n,'r-s')
grid on
xlabel('Frequency')
legend('Original Signal Spectrum','Sampled Signal Spectrum')
title('Frequency Domain of Original and Sampled Signal before LPF')

% Design a LPF for Reconstruction
f_cutoff = 4000;
f_stop = 6000;
lpFilt = designfilt('lowpassfir','PassbandFrequency',f_cutoff,'StopbandFrequency',f_stop,'Samplerate',Fs);
fvtool(lpFilt)
m1_rec = N*filter(lpFilt,s_out);

% we showed time domain
figure(4)
plot(t,m1,'b')
hold on
plot(t,m1_rec,'r')
grid on
legend('Original Message','Recovered Message')
title('Time Domain of Original and Recovered Signal')
xlabel('Time')
xlim([1.5 1.52])
sound(m1_rec,Fs)


% we showed freq. domain
fre_m1r = fftshift(fft(m1_rec,n));
figure(5)
stem(f,abs(fre_m1)/n,'b-s')
hold on
stem(f,abs(fre_m1r)/n,'r-o')
grid on
legend('Original Message Spectrum','Recovered Message Spectrum')
title('Frequency Domain of Original and Recovered Signal')
xlabel('Frequency')

%% Repeat the above steps for the new sampling rate 4 kHz. (1e)
% (b)
fs1 = 4000; % we sample the(approximately) cont time signal at a rate fs. every 5 samples for 8kHz
ts1 = 1/fs1;  % Sampling period
N1 = ts1/td; % N should be an integer N = 10
s_out1 = downsample(m1,N1);
s_out1 = upsample(s_out1,N1); 

% let's plot
figure(6)
plot(t,m1)
hold on
stem(t,s_out1,'r')
grid on
xlabel('Time')
legend('Original Signal','Sampled Signal','Location','SouthWest')
title('Original and Sampled Signal for 4 kHz')
xlim([2 2.001]) % very short time
% (c)
sa = dsp.SpectrumAnalyzer('SampleRate',Fs,'PlotAsTwoSidedSpectrum',true, ...
    'NumInputPorts',2,'ChannelNames', ...
    {'Message','Sampled Signal'});
sa(m1,s_out1);
release(sa);
% (d)
% take the Fourier Transform of m1 an s_out
f = (-(n-1)/2:(n-1)/2)*(Fs/n); % Generate the discrete frequency vector x axis 
fre_m1 = fftshift(fft(m1,n)); % Computes the FT of original (approximately) cont. time signal
fre_s_out1 = fftshift(fft(s_out1,n)); % Computes the FT of sampled disc. time signal
% let's plot
figure(7)
stem(f,(fre_m1)/n,'b-s')
hold on 
stem(f,(fre_s_out1)/n,'r-s')
grid on
xlabel('Frequency')
legend('Original Signal Spectrum','Sampled Signal Spectrum')
title('Frequency Domain of Original and Sampled Signal before LPF for 4 kHz')

% Design a LPF for Reconstruction
f_cutoff1 = 2000;
f_stop1 = 4000;
lpFilt1 = designfilt('lowpassfir','PassbandFrequency',f_cutoff1,'StopbandFrequency',f_stop1,'Samplerate',Fs);
fvtool(lpFilt1)
m1_rec1 = N1*filter(lpFilt1,s_out1);

% we showed time domain
figure(8)
plot(t,m1,'b')
hold on
plot(t,m1_rec1,'r')
grid on
legend('Original Message','Recovered Message')
title('Time Domain of Original and Recovered Signal for 4 kHz')
xlabel('Time')
xlim([1.5 1.52])
sound(m1_rec1,Fs)

% we showed freq. domain
fre_m1r = fftshift(fft(m1_rec1,n));
figure(9)
stem(f,abs(fre_m1)/n,'b-s')
hold on
stem(f,abs(fre_m1r)/n,'r-o')
grid on
legend('Original Message Spectrum','Recovered Message Spectrum')
title('Frequency Domain of Original and Recovered Signal for 4 kHz')
xlabel('Frequency')