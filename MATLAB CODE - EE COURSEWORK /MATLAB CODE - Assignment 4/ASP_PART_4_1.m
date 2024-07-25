% 4

% Define filter coefficients
b = [1 2 3 2 1];
a = 1;

% Generate WGN sequence x with 1000 samples
x = randn(1, 1000);

% Filter the WGN sequence using the filter function
y = filter(b, a, x);

% Normalize the output y so that its variance is unity
y = y/std(y);

% Generate noise sequence eta with standard deviation sigma = 0.1
eta = 0.1*randn(1, 1000);

% Add noise sequence to the filter output to obtain z
z = y + eta;

% Compute the SNR in dB
SNR_dB = 10*log10(sum(y.^2)/sum(eta.^2))


%%
% 4.1.1 


Rx = xcorr(x, 'biased');
Rx = Rx(ceil(length(Rx)/2):end);


Pzx = xcorr(z, x, 'biased');
Pzx = Pzx(ceil(length(Pzx)/2):end);

% Compute Rxx and pzx
Rxx = toeplitz(Rx);
pzx = Pzx';


w_opt = inv(Rxx)*pzx;


disp(w_opt);


disp(b);


%%
%4.1.1 - without normalization
clear all;
% Define filter coefficients
b = [1 2 3 2 1];
a = 1;


x = randn(1, 1000);


y = filter(b, a, x);


eta = 3.16*randn(1, 1000);


z = y + eta;


SNR_dB = 10*log10(sum(y.^2)/sum(eta.^2))


Rx = xcorr(x, 'biased');
Rx = Rx(ceil(length(Rx)/2):end);


Pzx = xcorr(z, x, 'biased');
Pzx = Pzx(ceil(length(Pzx)/2):end);

% Compute Rxx and pzx
Rxx = toeplitz(Rx);
pzx = Pzx';

% Compute optimal Wiener filter coefficients
w_opt = inv(Rxx)*pzx;

%%

clear all;
% Define filter coefficients
b = [1 2 3 2 1];
a = 1;

% Generate WGN sequence x with 1000 samples
x = randn(1, 1000);

% Filter the WGN sequence using the filter function
y = filter(b, a, x);

% Normalize the output y so that its variance is unity


% Generate noise sequence eta with standard deviation sigma = 0.1
eta = 0.01*randn(1, 1000);

% Add noise sequence to the filter output to obtain z
z = y + eta;

% Compute the SNR in dB
SNR_dB = 10*log10(sum(y.^2)/sum(eta.^2))

% Estimate autocorrelation function of x
Rx = xcorr(x, 'biased');
Rx = Rx(ceil(length(Rx)/2):end);

% Estimate cross-correlation function between x and z
Pzx = xcorr(z, x, 'biased');
Pzx = Pzx(ceil(length(Pzx)/2):end);

% Compute Rxx and pzx
Rxx = toeplitz(Rx);
pzx = Pzx';

% Compute optimal Wiener filter coefficients
w_opt = inv(Rxx)*pzx;

% Display optimal Wiener filter coefficients
disp(w_opt);

% Display unknown system filter coefficients
disp(b);