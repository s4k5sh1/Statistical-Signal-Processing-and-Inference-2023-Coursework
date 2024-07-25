 % 4. 5. 1 %chnage audio file to get plots for different sounds 


clear all;
fs = 44100; % Sampling frequency
N = 1000; % Number of samples
sound_files = {'sound_e.wav', 'sound_a.wav', 'sound_s.wav', 'sound_t.wav','sound_x.wav'};
sounds = {'e', 'a', 's', 't','x'};

% Different adaptation gains and predictor orders to investigate
mu_values = [0.0001, 0.001, 0.01, 0.1, 0.5, 1, 2, 5];
predictor_orders = [1, 2, 5, 10, 15, 20, 25, 30, 40, 50];

for s = 1:length(sounds)
    % Read the audio file
    [x, ~] = audioread(sound_files{s});
    x = x(1:N) / max(abs(x(1:N)));

    % Initialize MSE matrix
    mse_matrix = zeros(length(mu_values), length(predictor_orders));

    % Loop through the adaptation gains and predictor orders
    for i = 1:length(mu_values)
        mu = mu_values(i);
        for j = 1:length(predictor_orders)
            order = predictor_orders(j);

            % Apply the LMS algorithm
            [yhat, e, w] = lms(x(1:end-1), x(2:end), mu, order);

            % Compute the mean squared error (MSE)
            mse_matrix(i, j) = mean(e.^2);
        end
    end

    % Plot the heatmap
    figure;
    heatmap(predictor_orders, mu_values, mse_matrix, 'ColorScaling', 'log', 'CellLabelColor', 'none');
    title(['MSE for Sound "', sounds{s}, '"']);
    xlabel('Predictor Order');
    ylabel('Adaptation Gain');
    
    % Find the best adaptation gain and predictor order
    [min_mse, min_idx] = min(mse_matrix(:));
    [best_mu_idx, best_order_idx] = ind2sub(size(mse_matrix), min_idx);
    best_mu = mu_values(best_mu_idx);
    best_order = predictor_orders(best_order_idx);

end



 

%%
%4.5.3
 
clear all;
% Parameters
fs = 16000; % Sampling frequency
N = 5000; % Number of samples (increased)
sounds = {'e', 'a', 's', 't', 'x'};

% Different adaptation gains and predictor orders to investigate
mu_values = logspace(-4, 0, 10); % Logarithmically spaced adaptation gains
predictor_orders = [1, 2, 4, 8, 16, 32,64]; % Larger range of predictor orders

% Initialize prediction gains array
prediction_gains = zeros(1, length(sounds));

% Loop through each sound
for i = 1:length(sounds)
    % Read the audio file for the current sound
    [x, Fs] = audioread(['sound_', sounds{i}, '.wav']);
    x = x(1:N);

    x = x(1:N) / max(abs(x(1:N)));

    % Initialize MSE matrix
    mse_matrix = zeros(length(mu_values), length(predictor_orders));

    % Loop through the adaptation gains and predictor orders
    for j = 1:length(mu_values)
        mu = mu_values(j);
        for k = 1:length(predictor_orders)
            order = predictor_orders(k);

            % Apply the LMS algorithm
            [yhat, e, w] = lms(x(1:end-1), x(2:end), mu, order);

            % Compute the mean squared error (MSE)
            mse_matrix(j, k) = mean(e.^2);
        end
    end

    % Find the best adaptation gain and predictor order
    [min_mse, min_idx] = min(mse_matrix(:));
    [best_mu_idx, best_order_idx] = ind2sub(size(mse_matrix), min_idx);
    best_mu = mu_values(best_mu_idx)
    best_order = predictor_orders(best_order_idx)

    % Apply the LMS algorithm with the best adaptation gain and predictor order
    [yhat, e, w] = lms(x(1:end-1), x(2:end), best_mu, best_order);

    % Calculate the prediction gain
    input_variance = var(x);
    error_variance = var(e);
    prediction_gains(i) = 10 * log10(input_variance / error_variance);
end

% Display the prediction gains
disp(['Prediction gains: ', num2str(prediction_gains)]);


%%
function [yhat, e, w] = lms(x, z, mu, Nw)

% x: input signal
% z: desired signal
% mu: adaptation gain
% Nw: filter order

% Initialization
N = length(x);          % Length of input signal
w = zeros(Nw+1, N);     % Matrix containing the evolution of the adaptive weights in time
yhat = zeros(N, 1);     % Estimated output signal
e = zeros(N, 1);        % Error signal

% LMS algorithm
for n = 1:N
    xvec = [flipud(x(max(n-Nw,1):n)); zeros(Nw-n+1, 1)]; % Construct input vector
    yhat(n) = w(:,n)' * xvec;  % Compute estimated output
    e(n) = z(n) - yhat(n);     % Compute error
    w(:,n+1) = w(:,n) + mu * e(n) * xvec; % Update filter coefficients
end

end