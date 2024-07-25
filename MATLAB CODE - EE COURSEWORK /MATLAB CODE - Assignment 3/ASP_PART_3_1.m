% pgm - 3 
clear all;
N = 128;
x = randn(N, 1);
[PbX, freq] = pgm_norm(x);
subplot(3,1,1)
plot(freq, PbX);
title('Periodogram of WGN for N = 128');
xlabel('Normalized frequency (× 2π rad/sample)');
ylabel('Power spectral density');

N = 256;
x = randn(N, 1);
[PbX, freq] = pgm_norm(x);
subplot(3,1,2)
plot(freq, PbX);
title('Periodogram of WGN for N = 256');
xlabel('Normalized frequency (× 2π rad/sample)');
ylabel('Power spectral density');

N = 512;
x = randn(N, 1);
[PbX, freq] = pgm_norm(x);
subplot(3,1,3)
plot(freq, PbX);
title('Periodogram of WGN for N = 512');
xlabel('Normalized frequency (× 2π rad/sample)');
ylabel('Power spectral density');

%%
% 3.1.1

clear all;
N = 128;
x = randn(N, 1);
[PbX, freq] = pgm_norm(x);


b = 0.2 * [1 1 1 1 1];
PbX_filtered = filter(b, 1, PbX);


f = linspace(0, 1, N);
subplot(3,1,1)
plot(f, PbX, 'b', f, PbX_filtered, 'r');
title('Periodogram of WGN for N = 128');
xlabel('Normalized frequency (× 2π rad/sample)');
ylabel('Power spectral density');
legend('Original PSD Estimate', 'Smoothed PSD Estimate');

N = 256;
x = randn(N, 1);
[PbX, freq] = pgm_norm(x);

b = 0.2 * [1 1 1 1 1];
PbX_filtered = filter(b, 1, PbX);


f = linspace(0, 1, N);
subplot(3,1,2)
plot(f, PbX, 'b', f, PbX_filtered, 'r');
title('Periodogram of WGN for N = 256');
xlabel('Normalized frequency (× 2π rad/sample)');
ylabel('Power spectral density');
legend('Original PSD Estimate', 'Smoothed PSD Estimate');

N = 512;
x = randn(N, 1);
[PbX, freq] = pgm_norm(x);

b = 0.2 * [1 1 1 1 1];
PbX_filtered = filter(b, 1, PbX);


f = linspace(0, 1, N);
subplot(3,1,3)
plot(f, PbX, 'b', f, PbX_filtered, 'r');
title('Periodogram of WGN for N = 512');
xlabel('Normalized frequency (× 2π rad/sample)');
ylabel('Power spectral density');
legend('Original PSD Estimate', 'Smoothed PSD Estimate');


%%
% 3.1.2
clear all;

x = randn(1024, 1);


segments = reshape(x, [128, 8]);


psd_estimates = zeros(128, 8);
for i = 1:8
    [psd_estimates(:, i), freq] = pgm_norm(segments(:, i));
end


figure;

for i = 1:8
    subplot(4,2,i);
    plot(freq, psd_estimates(:, i));
    xlabel('Normalized frequency (× 2π rad/sample) ');
    ylabel('PSD');
    title(sprintf('Segment %d',i));
end


%%
clear all;

x = randn(1024, 1);

% Divide x into 8 non-overlapping 128-sample segments
segments = reshape(x, 128, 8);


psd_segments = zeros(128, 8);
for i = 1:8
    [psd_segments(:,i), freq] = pgm_norm(segments(:,i));
end

%
avg_psd = mean(psd_segments, 2);


plot(freq, avg_psd, 'b');
xlabel('Normalized frequency (× 2π rad/sample)');
ylabel('Power spectral density');
title('Averaged periodogram');


%%
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