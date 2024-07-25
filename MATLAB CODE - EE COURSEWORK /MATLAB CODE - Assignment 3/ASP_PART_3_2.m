% 3.2.1 
clear all;
[h,w]=freqz([1],[1 0.9],512);
figure;
plot(w/(2*pi),abs(h).^2)
xlabel('Normalized frequency (× 2π rad/sample)');
ylabel('Power spectral density');
title('Exact PSD of AR(1) Process')


% 3.2.2

clear all;
% Generate WGN sequence
x = randn(1, 1064);
% Filter the signal
b = 1;
a = [1, 0.9];
y = filter(b, a, x);

[h,w]=freqz([1],[1 0.9],512);
figure;
plot(w/(2*pi),abs(h).^2,'b',LineWidth=2)
hold on;
[PbX, freq] = pgm_norm(y);
plot(freq, PbX, 'r');
hold off;
xlabel('Normalized frequency (× 2π rad/sample)');
ylabel('Power spectral density');
title('PSD of AR(1) Process')
xlim([0 0.5])
ylim([0 150])
legend(["Ideal PSD","Periodogram Estimated PSD"])

%3.2.3

clear all;
% Generate WGN sequence
x = randn(1, 1064);
% Filter the signal
b = 1;
a = [1, 0.9];
y = filter(b, a, x);

[h,w]=freqz([1],[1 0.9],512);
figure;
plot(w/(2*pi),abs(h).^2,'b',LineWidth=2)
hold on;
[PbX, freq] = pgm_norm(y);
plot(freq, PbX, 'r');
hold off;
xlabel('Normalized frequency (× 2π rad/sample)');
ylabel('Power spectral density');
title('PSD of AR(1) Process')
xlim([0.4 0.5])
legend(["Ideal PSD","Periodogram Estimated PSD"])



%3.2.4


clear all;
% Generate WGN sequence
x = randn(1, 1064);
% Filter the signal
b = 1;
a = [1, 0.9];
y = filter(b, a, x);

N = length(y);
[r, lags] = xcorr(y, 'biased');
r = r(N:end);
a1 = -r(2)/r(1);
sigma2 = r(1) + a1*r(2);
[h, w] = freqz(sigma2, [1, a1], 512);
Pby = abs(h).^2;

[h,w]=freqz([1],[1 0.9],512);
figure;
plot(w/(2*pi),abs(h).^2,'b',LineWidth=2)
hold on;
[PbX, freq] = pgm_norm(y);
plot(freq, PbX, 'r',w/(2*pi), Pby, 'k',LineWidth=2);
hold off;
xlabel('Normalized frequency (× 2π rad/sample)');
ylabel('Power spectral density');
title('Model Based PSD Estimate')
xlim([0.4 0.5]);
legend(["Ideal PSD","Periodogram Estimated PSD","Model Based Estimated PSD"])


% 3.2.5

clear all;
load sunspot.dat
y = sunspot(:,2);

% periodoogram 
[Pby, freq] = pgm_norm(y);
plot(freq, Pby);
xlabel('Normalized Frequency (× 2π rad/sample)');
ylabel('PSD');
title('Periodogram of Original Data')
xlim([0 0.2]);

% periodogram and model based psd for mean centered data 
y_centered = y - mean(y);
[Pby_centered, freq] = pgm_norm(y_centered);
order=10;
subplot(1,4,2)
plot(freq, Pby_centered);
xlabel('Normalized Frequency (× 2π rad/sample)');
ylabel('PSD');
xlim([0 0.5]);
subplot(1,4,3)
[a_centered, sigma_centered] = aryule(y_centered, order);
[Pby_model_centered, w] = freqz(sigma_centered, a_centered, 512);
plot(w/(2*pi), abs(Pby_model_centered).^2);
xlabel('Normalized Frequency (× 2π rad/sample)');
ylabel('PSD');
xlim([0 0.2]);


% model based psd for different model orders
orders = [1, 5, 10, 20, 50]; 
subplot(1,4,4)
hold on;
plot(freq, Pby);
for i = 1:length(orders)
    order = orders(i);
    [a, sigma] = aryule(y, order);
    b = 1;
    [Pby_model, w] = freqz(sigma, a, 512);
    plot(w/(2*pi), abs(Pby_model).^2);
end
hold off;
xlabel('Normalized Frequency');
ylabel('PSD');
xlim([0 0.5]);
legend('Periodogram', 'Order 1', 'Order 5', 'Order 10', 'Order 20', 'Order 50');

function [PbX, freq] = pgm_norm(x)
% Calculates the periodogram of a sequence x using the given equation
% with a normalized frequency axis going from 0 to 1
%   PbX: the periodogram of x, a sequence of length N
%   freq: the normalized frequency axis, a sequence of length N

N = length(x);
PbX = zeros(N, 1);
freq = (0:N-1)/N; % normalized frequency axis
for f = 1:N
    for n = 1:N
        PbX(f) = PbX(f) + x(n) * exp(-1i*2*pi*(f-1)*(n-1)/N);
    end
    PbX(f) = abs(PbX(f))^2 / N;
end
end 

