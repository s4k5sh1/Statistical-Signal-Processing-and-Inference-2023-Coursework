% %3.4.1
clear all;
% Generate a random London landline number
number = sprintf('020 %04d %04d', randi([0 9], [1 2 4]));

disp(number);

f = [697 1209;
     697 1336;
     697 1477;
     770 1209;
     770 1336;
     770 1477;
     852 1209;
     852 1336;
     852 1477;
     941 1209;
     941 1336;
     941 1477];

% Define the DTMF symbols
symbols = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '*', '0', '#'];

% Compute the DTMF signals
duration = 0.25;
fs = 32768;
t = 0:1/fs:duration-1/fs;
dtmf_signals = zeros(length(symbols), length(t));
for i = 1:length(symbols)
    f1 = f(i, 1);
    f2 = f(i, 2);
    signal = sin(2*pi*f1*t) + sin(2*pi*f2*t);
    dtmf_signals(i, :) = signal;
end

% Concatenate the signals for digit 0 and 2 with an idle time in between
signal_0 = dtmf_signals(11, :);
idle_time = zeros(1, length(t));
signal_2 = dtmf_signals(2, :);
dtmf_sequence = [signal_0 idle_time signal_2];

% Compute the time vector
duration = 0.75;
t = 0:1/fs:duration-1/fs;

% Plot the DTMF sequence
figure
plot(t, dtmf_sequence)
xlabel('Time (s)')
ylabel('Amplitude')
title('DTMF Sequence')
xlim([0 0.75]);

% 3.4.2
f = [697 1209;
     697 1336;
     697 1477;
     770 1209;
     770 1336;
     770 1477;
     852 1209;
     852 1336;
     852 1477;
     941 1209;
     941 1336;
     941 1477];

% Define the DTMF symbols
symbols = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '*', '0', '#'];

% Compute the DTMF signals
duration = 0.25;
fs = 32768;
t = 0:1/fs:duration-1/fs;
dtmf_signals = zeros(length(symbols), length(t));
for i = 1:length(symbols)
    f1 = f(i, 1);
    f2 = f(i, 2);
    signal = sin(2*pi*f1*t) + sin(2*pi*f2*t);
    dtmf_signals(i, :) = signal;
end

% Concatenate the signals for digit 0 and 2 with an idle time in between
signal_0 = dtmf_signals(11, :);
idle_time = zeros(1, length(t));
signal_2 = dtmf_signals(2, :);
dtmf_sequence = [signal_0 idle_time signal_2];

% Compute the time vector
duration = 0.75;
t = 0:1/fs:duration-1/fs;

% Plot the DTMF sequence
figure
plot(t, dtmf_sequence)
xlabel('Time (s)')
ylabel('Amplitude')
title('DTMF Sequence')
xlim([0 0.75]);


% Define spectrogram parameters
nfft = 8192;
window = hann(nfft);

% Compute spectrogram
[~, f, t, P] = spectrogram(dtmf_sequence, window, 0, nfft, fs, 'yaxis');

% Scale frequency vector to kHz
f = f/1000;

% Plot spectrogram
figure
imagesc(t, f, 10*log10(P))
set(gca, 'YDir', 'normal')
xlabel('Time (s)')
ylabel('Frequency (kHz)')
title('DTMF Spectrogram')
ylim([0 1.6])
colorbar







 
% 3.4.4 - variance - 5
% Generate random phone number
phone_num = sprintf('020%04d%04d', randi(9999), randi(9999));
disp(phone_num)

% Convert phone number to DTMF sequence
dtmf_sequence = [];
for i = 1:length(phone_num)
    digit = phone_num(i);
    idx = find(symbols == digit);
    if ~isempty(idx)
        signal = dtmf_signals(idx, :);
        dtmf_sequence = [dtmf_sequence signal];
        % Add idle time with noise between each digit
        if i < length(phone_num)
            idle_time = zeros(1, round(0.25*fs));
            noise = sqrt(variance)*randn(1, length(idle_time));
            noisy_idle_time = idle_time + noise;
            dtmf_sequence = [dtmf_sequence noisy_idle_time];
        end
    end
end

% Add white Gaussian noise to the sequence
variance = 5;
noise = sqrt(variance)*randn(1, length(dtmf_sequence));
noisy_dtmf_sequence = dtmf_sequence + noise;

% Compute time vector
duration = length(noisy_dtmf_sequence)/fs;
t = 0:1/fs:duration-1/fs;

% Plot full time domain signal
figure
subplot(2,1,1)
plot(t, noisy_dtmf_sequence)
xlabel('Time (s)')
ylabel('Amplitude')
title('Noisy DTMF Signal')
xlim([0 5.25])

% Compute spectrogram
nfft = 4096;
window = hann(nfft);
[S, f, t] = spectrogram(noisy_dtmf_sequence, window, 0, nfft, fs, 'yaxis');

% Scale frequency vector to kHz
f = f/1000;

% Plot spectrogram
subplot(2,1,2)
surf(t, f, 10*log10(abs(S)), 'EdgeColor', 'none')
view(2)
xlabel('Time (s)')
ylabel('Frequency (kHz)')
title('Noisy DTMF Spectrogram')
ylim([0 1.6])
colorbar


% 3.4.4 - variance - 50
% Generate random phone number
phone_num = sprintf('020%04d%04d', randi(9999), randi(9999));
disp(phone_num)

% Convert phone number to DTMF sequence
dtmf_sequence = [];
for i = 1:length(phone_num)
    digit = phone_num(i);
    idx = find(symbols == digit);
    if ~isempty(idx)
        signal = dtmf_signals(idx, :);
        dtmf_sequence = [dtmf_sequence signal];
        % Add idle time with noise between each digit
        if i < length(phone_num)
            idle_time = zeros(1, round(0.25*fs));
            noise = sqrt(variance)*randn(1, length(idle_time));
            noisy_idle_time = idle_time + noise;
            dtmf_sequence = [dtmf_sequence noisy_idle_time];
        end
    end
end

% Add white Gaussian noise to the sequence
variance = 50;
noise = sqrt(variance)*randn(1, length(dtmf_sequence));
noisy_dtmf_sequence = dtmf_sequence + noise;

% Compute time vector
duration = length(noisy_dtmf_sequence)/fs;
t = 0:1/fs:duration-1/fs;

% Plot full time domain signal
figure
subplot(2,1,1)
plot(t, noisy_dtmf_sequence)
xlabel('Time (s)')
ylabel('Amplitude')
title('Noisy DTMF Signal')
xlim([0 5.25])

% Compute spectrogram
nfft = 4096;
window = hann(nfft);
[S, f, t] = spectrogram(noisy_dtmf_sequence, window, 0, nfft, fs, 'yaxis');

% Scale frequency vector to kHz
f = f/1000;

% Plot spectrogram
subplot(2,1,2)
surf(t, f, 10*log10(abs(S)), 'EdgeColor', 'none')
view(2)
xlabel('Time (s)')
ylabel('Frequency (kHz)')
title('Noisy DTMF Spectrogram')
ylim([0 1.6])
colorbar


% 3.4.4 - variance - 500
% Generate random phone number
phone_num = sprintf('020%04d%04d', randi(9999), randi(9999));
disp(phone_num)

% Convert phone number to DTMF sequence
dtmf_sequence = [];
for i = 1:length(phone_num)
    digit = phone_num(i);
    idx = find(symbols == digit);
    if ~isempty(idx)
        signal = dtmf_signals(idx, :);
        dtmf_sequence = [dtmf_sequence signal];
        % Add idle time with noise between each digit
        if i < length(phone_num)
            idle_time = zeros(1, round(0.25*fs));
            noise = sqrt(variance)*randn(1, length(idle_time));
            noisy_idle_time = idle_time + noise;
            dtmf_sequence = [dtmf_sequence noisy_idle_time];
        end
    end
end

% Add white Gaussian noise to the sequence
variance = 500;
noise = sqrt(variance)*randn(1, length(dtmf_sequence));
noisy_dtmf_sequence = dtmf_sequence + noise;

% Compute time vector
duration = length(noisy_dtmf_sequence)/fs;
t = 0:1/fs:duration-1/fs;

% Plot full time domain signal
figure
subplot(2,1,1)
plot(t, noisy_dtmf_sequence)
xlabel('Time (s)')
ylabel('Amplitude')
title('Noisy DTMF Signal')
xlim([0 5.25])

% Compute spectrogram
nfft = 4096;
window = hann(nfft);
[S, f, t] = spectrogram(noisy_dtmf_sequence, window, 0, nfft, fs, 'yaxis');

% Scale frequency vector to kHz
f = f/1000;

% Plot spectrogram
subplot(2,1,2)
surf(t, f, 10*log10(abs(S)), 'EdgeColor', 'none')
view(2)
xlabel('Time (s)')
ylabel('Frequency (kHz)')
title('Noisy DTMF Spectrogram')
ylim([0 1.6])
colorbar