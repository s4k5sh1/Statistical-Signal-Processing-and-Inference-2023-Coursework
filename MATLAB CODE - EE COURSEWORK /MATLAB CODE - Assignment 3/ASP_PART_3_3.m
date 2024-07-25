% 3.3.3 
% Load the sunspot data
load sunspot.dat
sunspot_data = sunspot(:, 2);

% Normalize the data
sunspot_data = (sunspot_data - mean(sunspot_data)) / std(sunspot_data);

% Initialize the parameters
max_order = 10;
n = length(sunspot_data);

% Preallocate arrays
AR_coeffs = cell(max_order, 1);
errors = zeros(max_order, 1);

% Loop over model orders
for p = 1:max_order
    % Calculate autocorrelation
    rxx = xcorr(sunspot_data, p, 'biased');
    rxx = rxx(p+1:end);
    
    % Prepare the H matrix and x vector
    H = toeplitz(rxx(1:end-1));
    x = rxx(2:end);
    
    % Calculate the LSE coefficients
    a = H \ x;
    
    % Store the coefficients and errors
    AR_coeffs{p} = a;
    errors(p) = norm(x - H * a) ^ 2;
end

% Display the results
disp(AR_coeffs);
disp(errors);

figure;
hold on;

for p = 1:max_order
    subplot(5, 2, p);
    bar(AR_coeffs{p});
    title(['AR(' num2str(p) ') coefficients']);
    xlabel('Model Order');
    ylabel('Value');
    ylim([-1, 1]);
end

hold off;

%%
%3.3.4

figure;
plot(1:max_order, errors, '-o');
title('Approximation Error vs. Model Order');
xlabel('Model Order (p)');
ylabel('Squared Error');

%%
% 3.3.5

figure;
for p = 1:max_order
    % Extract AR(p) coefficients
    ap = AR_coeffs{p};
    
    % Calculate the AR(p) model
    ARp_data = filter([1 -ap'], 1, sunspot_data);
    
  
    
    % Compute the power spectral density (PSD) using the periodogram function
    %[psd_ARp, freq] = periodogram(ARp_data, [], [], 1);
     [PbX, freq] = pgm_norm(ARp_data);
    % Plot the power spectra for the current order
    subplot(5, 2, p);
    plot(freq, PbX);
    title(['PSD of AR(' num2str(p) ') Model - Sunspot Time Series']);
    xlabel('Normalized Frequency');
    ylabel('PSD');
    
end




%%
% 3.3.6

% Define data lengths and model order
load sunspot.dat
sunspot_data = sunspot(:, 2);

% Normalize the data
sunspot_data = (sunspot_data - mean(sunspot_data)) / std(sunspot_data);

data_lengths = [10: 5: 250];
model_order = 2;

% Initialize error array
errors_N = zeros(1, length(data_lengths));

% Loop over data lengths
for idx = 1:length(data_lengths)
    N = data_lengths(idx);
    
    % Take the first N samples of the sunspot_data
    sunspot_data_N = sunspot_data(1:N);
    
    % Calculate autocorrelation
    rxx = xcorr(sunspot_data_N, model_order, 'unbiased');
    rxx = rxx(model_order+1:end);
    
    % Prepare the H matrix and x vector
    H = toeplitz(rxx(1:end-1));
    x = rxx(2:end);
    
    % Calculate the LSE coefficients
    a = H \ x;
    
    % Store the approximation error
    errors_N(idx) = norm(x - H * a) ^ 2;
end

% Plot the approximation error
figure;
bar(data_lengths, errors_N);
ylabel('Squared Error')
xlabel('Data Length (N)')
title('Approxiamtion error for different data length N')






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