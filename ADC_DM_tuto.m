%%%%%%%%%%%Sinus Generation%%%%%%%%%%%%%%
t_sim = 0:0.01:10 ; % Time is inherently sampled
y = cos(2*pi*t_sim);  % Cosinus is also sampled
plot(t_sim,y)



%%%%%%%%%%%Sinus Generation -discrete%%%%%%%%%%%%%%

t_sim = 0:0.1:3;
y  = cos(2*pi*t_sim);
% Prepare figure with two plots
subplot(211) % Use stem to display the sampled sequence
stem(t_sim,y,'linewidth',2)
xlabel('time')
ylabel('y')

% Enable second plot
subplot(212) % Use plot to display the sampled sequence
plot(t_sim,y,':o')
xlabel('time')
ylabel('y')

%%%%%%%%%%% Quantization Code%%%%%%%%%%%%%%
t_sim = 0:0.01:2;
Vref = 1.5; Nbits = 4;
x = 1.35*cos(2*pi*1*t_sim);
quantizedInput = floor((x+Vref)*2^(Nbits-1)/Vref); % Quantizing the sampled data
quantizedInput(quantizedInput<0)  = 0; % Clipping Down 
quantizedInput(quantizedInput>2^Nbits-1) = 2^Nbits-1; % Clipping Up
DigOutput = (quantizedInput-2^(Nbits-1))/2^(Nbits-1)*Vref+Vref/2^Nbits;
stem(t_sim,DigOutput) ; xlabel('Temps (s)') ; ylabel('Sortie Quantifiee')

%%%%%%%%%%% Up Sampling %%%%%%%%%%%%%%
fs = 10;
tstop = 1.75;
t = 0:1/fs:tstop;
f = 1;
y = cos(2*pi*f*t);
USR  = 4; % upsampling ratio
fs_up  = fs*USR;
y_up = zeros(1,(length(y)-1)*USR+1);
y_up(1:USR:end) = y;
t_up = 0:1/fs_up:t(end);

%%%%%%%%%%% Up Sampling %%%%%%%%%%%%%%
fs = 10;
tstop = 1.75;
t = 0:1/fs:tstop;
f = 1;
y = cos(2*pi*f*t);
DSR = 4;
fs_down = fs/DSR;
y_down = y(1:DSR:end);
t_down = 0:1/fs_down:tstop;




%%%%%%%%%%% Random signal spectrum %%%%%%%%%%%%%%
Fs = 1;
x = rand(51,1) + 1i*rand(51,1); % Complex signal
Xpsd = abs(fft(x)).^2; % Note the dot ! Note the square !
Nx  = length(Xpsd); % length of the FFT, also length of x
bin_freq_val = [0:Nx-1];
subplot(2,1,1); stem(0:Nx-1,real(x))
xlabel('Time index');
 ylabel('Magnitude')
subplot(2,1,2); plot(bin_freq_val,10*log10(Xpsd)) % Note the 10xlog10 !
xlabel('Frequency bin'); ylabel('Power spectral density (dB)')

figure()
bin_freq_val_shift = -(Nx-1)/2 : (Nx-1)/2;
freq_val_shift = bin_freq_val_shift/Nx*Fs;
plot(freq_val_shift,fftshift(10*log10(Xpsd))) ; xlim(0.5*[-1 1])
xlabel('Frequency (Hz)'); ylabel('Power spectral density (dB)')

%%%%%%%%%%%%%%%%%%%%%%% Power and error Calculation basics%%%%%%%%%%%%%%%%%
clear;
close all;

Nsim    = 2^14;                         % Number of points of the simulation
Fs      = 30e6;                         % Sampling Frequency
Ts      = 1/Fs;                         % Sampling Period
t_sim   = 0:Ts:(Nsim-1)*Ts;             % Time Vector
t_sim   = reshape(t_sim,Nsim,1);        % Reshape Row vector to Column vector

fsig_or   = 1e6;                        % Input frequency (Hz)
sig_bin   = round(fsig_or/Fs*Nsim);     % Input bin of the signal frequency
fsig      = sig_bin*Fs/Nsim;            % re-adjusting in a signal bin

Amp   = 1;                              % Signal Amplitude
x     = Amp*sin(2*pi*fsig*t_sim);       % Input signal generation

sigma_noise   = 0.01;                   % Sigma of the noise (Ecart type)
x_noi         = x+sigma_noise*randn(size(t_sim));  % Adding noise to the signal
error         = x_noi - x;              % Calculating the error

PS_theo   = Amp.^2/2;                   % Theoretical Average Signal power
PN_theo   = sigma_noise.^2;             % Theoretical error or noise power
SNR_theo  = 10*log10(PS_theo/PN_theo);  % Theoretical SNR
disp(['The theoretical SNR is ', num2str(SNR_theo), ' dB'])

PS_time   = mean(x.^2);                 % Empirical time domain Signal power 
PN_time   = mean(error.^2);             % Empirical time domain error or noise power 
SNR_time  = 10*log10(PS_time/PN_time);  % Empirical time domain SNR 
disp(['The time domain SNR is ', num2str(SNR_time), ' dB'])

Nx_noi    = length(x_noi);
win       = blackman(Nx_noi,'periodic');   % Calculating the window
x_noiPSD  = abs(fft(x_noi(:).*win(:))).^2; % Calculating the PSD other windowed signal
x_noiPSD  = x_noiPSD/Nx_noi;               % Normalizing with respect to the length

figure()
plot(10*log10(x_noiPSD))
xlabel('bin index')
ylabel('PSD(dB/bin)') 
title('Plot of the raw FFT');

figure()
subplot(2,1,1)
plot(10*log10(x_noiPSD(1:Nsim/2)))
xlabel('bin index')
ylabel('PSD(dB/bin)') 
title('Plot of the FFT (positive frequencies only)');

subplot(2,1,2)
f_step  = Fs/Nsim;
f       = 0:f_step:(Nsim/2-1)*f_step;
plot(f/1e6,10*log10(x_noiPSD(1:Nsim/2)))
xlabel('Frequency (MHz)')
ylabel('PSD(dB/bin)') 


err_bin   = [1 , Nsim/2];   % beggining and end of the useful band (here, Nyquist band)
nintsig   = 2;              % Number of points around the useful signal to integrate
PS_freq   = sum(x_noiPSD(sig_bin+1-nintsig:sig_bin+1+nintsig)); % Empirical frequency domain Signal power 
PN_freq   = sum(x_noiPSD(1:Nsim/2))-PS_freq;                    % Empirical frequency domain Noise or error power
SNR_freq  = 10*log10(PS_freq/PN_freq);
disp(['The frequency domain SNR is ', num2str(SNR_freq), ' dB'])


