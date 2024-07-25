clear all;
figure(6);
hold on
xlim([-2.5 2.5])
ylim([-1.5 1.5])
xlabel('a_1')
ylabel('a_2')
title('Stability Triangle')

for i = 1:1000
    a1 = unifrnd(-2.5, 2.5);
    a2 = unifrnd(-1.5, 1.5);
    roots_func = roots([1, -a1,-a2]);
    
    % Generate AR(2) process with coefficients a1 and a2
    x = zeros(1000,1);
    x(1) = randn;
    x(2) = randn;
    for n = 3:1000
        x(n) = a1*x(n-1) + a2*x(n-2) + randn;
    end
    
    % Check for stable
    if all(abs(roots_func) < 1)
        plot(a1, a2, '*b')
    else
        plot(a1, a2, '*r')
    end
end

%%
% 2.3.2

clear all;
data = load('sunspot.dat');
sunspot = data(:,2); %sunspot data from the second column

N = 5; % Data length
% Compute ACF of sunspot data up to lag N
[acf, lags] = xcorr(sunspot(1:N), 'unbiased');
figure(9);
% Plot ACF
subplot(3,1,1)
stem(lags, acf);
xlabel('Correlation lag');
ylabel('Correlation');
title('ACF for Sunspot Data with N=5');

N =20 ;
[acf, lags] = xcorr(sunspot(1:N), 'unbiased');
subplot(3,1,2)
stem(lags, acf);
xlabel('Correlation lag');
ylabel('Correlation');
title('ACF for Sunspot Data with N=20');

N =250 ;
[acf, lags] = xcorr(sunspot(1:N), 'unbiased');
subplot(3,1,3)
stem(lags, acf);
xlabel('Correlation lag');
ylabel('Correlation');
title('ACF for Sunspot Data with N=250');

%%
% 2.3.2 - normalised version for comparision 
clear all;
load sunspot.dat % Load sunspot data
sunspot_data = sunspot(:, 2); % Extract sunspot data from second column

N = 5; % Data length
sunspot_data = sunspot_data(1:N) - mean(sunspot_data(1:N)); % Zero-mean version of sunspot data up to N

% Compute ACF of zero-mean sunspot data up to lag N-1
[acf, lags] = xcorr(sunspot_data, N-1, 'coeff');

% Plot ACF
subplot(3,1,1)
stem(lags, acf);
xlabel('Lag');
ylabel('ACF');
title(sprintf('ACF for Zero-Mean Sunspot Data (N=%d)', N));

load sunspot.dat % Load sunspot data
sunspot_data = sunspot(:, 2); % Extract second column of data
N = 20; % Data length

% Compute zero-mean version of sunspot data
sunspot_zero_mean = sunspot_data(1:N) - mean(sunspot_data(1:N));

% Compute ACF of zero-mean sunspot data up to lag N
[acf, lags] = xcorr(sunspot_zero_mean(1:N), N-1, 'coeff');

% Plot ACF
subplot(3,1,2)
stem(lags, acf);
xlabel('Lag');
ylabel('ACF');
title(sprintf('ACF for Zero-Mean Sunspot Data (N=%d)', N));

load sunspot.dat % Load sunspot data
sunspot_data = sunspot(:, 2); % Extract second column of data
N = 250; % Data length

% Compute zero-mean version of sunspot data
sunspot_zero_mean = sunspot_data(1:N) - mean(sunspot_data(1:N));

% Compute ACF of zero-mean sunspot data up to lag N
[acf, lags] = xcorr(sunspot_zero_mean(1:N), N-1, 'coeff');

% Plot ACF
subplot(3,1,3)
stem(lags, acf);
xlabel('Lag');
ylabel('ACF');
title(sprintf('ACF for Zero-Mean Sunspot Data (N=%d)', N));


%%
% 2.3.3
clear all;
load sunspot.dat; % Load sunspot data
x = sunspot(:,2); % Extract data from second column
[arcoefs,E,K] = aryule(x,15);
pacf = -K;
stem(pacf)
xlabel('Lag')
ylabel('Partial ACF')
title('Partial Autocorrelation Sequence')
xlim([1 10])
hold on
x_std = (x - mean(x)) / std(x);
[arcoefs,E,K] = aryule(x_std,15);
pacf = -K;
stem(pacf)
xlabel('Lag')
ylabel('Partial ACF')
title('Partial Autocorrelation Sequence')
conf = sqrt(2)*erfinv(0.95)/sqrt(1000);
plot(xlim,[1 1]'*[-conf conf],'r')
hold off
legend(["Empirical PCF","Standardised PCF"])

%%
% 2.3.4

clear all;
load sunspot.dat; 
sunspot_data = sunspot(:,2); 

% Standardize the data to zero mean and unit variance
sunspot_std = (sunspot_data - mean(sunspot_data)) / std(sunspot_data);


max_p = 10; % Maximum model order
mdl = zeros(max_p, 1);
aic = zeros(max_p, 1);
aicc = zeros(max_p, 1);

for p = 1:max_p
    
    [a, e, k] = aryule(sunspot_std, p);
    

    n = length(sunspot_std);
    mdl(p) = log(e) + p*log(n)/n;
    aic(p) = log(e) + 2*p/n;
    aicc(p) = aic(p) + (2*p*(p+1))/(n-p-1);
end

% Plot the MDL, AIC, and AICc values for different model orders
figure;
plot(1:max_p, mdl, 'o-', 1:max_p, aic, 'o-', 1:max_p, aicc, 'o-');
legend('MDL', 'AIC', 'AICc');
title('MDL,AIC and AICc from p = 0 to p = 10')
xlabel('Model Order');
ylabel('Critical Value');




%%
% 2.3.5 - ar modelling sunspot data % change prediction horizon values to
% get graph 

clear all;
load sunspot.dat; % Load sunspot data
sunspot_series = sunspot(:, 2);

% Standardize the data to zero mean and unit variance
sunspot_series_standardized = (sunspot_series - mean(sunspot_series)) / std(sunspot_series);

% Define the prediction horizon
prediction_horizon = 10;

% AR model orders to test
model_orders = [1, 2, 10];

% Number of orders
num_orders = length(model_orders);

% Prepare matrices to store the predicted values
predicted_values_matrix = zeros(length(sunspot_series_standardized), num_orders);

% Loop through the model orders
for i = 1:num_orders
    order = model_orders(i);

    % Calculate AR model coefficients using the aryule function
    [AR_coeffs, ~] = aryule(sunspot_series_standardized, order);

    % Calculate the predictions
    predicted_values = filter(-AR_coeffs(2:end), 1, sunspot_series_standardized);
    predicted_values = [zeros(prediction_horizon, 1); predicted_values(1:end-prediction_horizon)];

    % Store the predicted values in the matrix
    predicted_values_matrix(:, i) = predicted_values;
end

% Plot the actual and predicted values in subplots
figure;
for i = 1:num_orders
    order = model_orders(i);

    subplot(1, num_orders, i);
    plot(sunspot_series_standardized, 'LineWidth', 1.5);
    hold on;
    plot(predicted_values_matrix(:, i), 'r');
    hold off;

    title(sprintf('AR(%d)', order));
    xlabel('Time');
    ylabel('Standardized Sunspot Number');
    xlim([0 100]);   
end
legend('Actual Standardized Sunspot Data','AR model for prediciton horizon m=2')