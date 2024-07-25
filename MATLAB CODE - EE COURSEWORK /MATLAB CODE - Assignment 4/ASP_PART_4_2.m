%4.2.1 

% Generate input signal x
N = 1000;
x = randn(N, 1);

% Define unknown system
b = [1, 2, 3, 2, 1];
a = 1;

% Filter input signal to get unknown system output y
y = filter(b, a, x);


y = y / std(y);


sigma = 0.1;
eta = sigma * randn(N, 1);


z = y + eta;

% Apply LMS algorithm
mu = 0.01;   % Adaptation gain
Nw = 4;      % Filter order
[yhat, e, w] = lms(x, z, mu, Nw);


figure;
subplot(4,1,1);
plot(y);
title('Unknown system output : y[n]');
xlabel('Time');
ylabel('Signal Amplitude');
subplot(4,1,2);
plot(z);
title('Desired signal : z[n]');
xlabel('Time');
ylabel('Signal Amplitude');
subplot(4,1,3);
plot(yhat);
title('LMS estimate');
xlabel('Time');
ylabel('Signal Amplitude');
subplot(4,1,4);
plot(yhat,'b',LineWidth=1);
title('LMS estimate and Desired signal : z[n] ');
xlabel('Time');
ylabel('Signal Amplitude');
hold on;
plot(z);
legend('LMS Estimate','Signal z[n]')
hold off;



% 4.2.2 
% chnage mu to get different plots 
% Generate input signal x
N = 1000;
x = randn(N, 1);

% Define unknown system
b = [1, 2, 3, 2, 1];
a = 1;

% Filter input signal to get unknown system output y
y = filter(b, a, x);

% Normalize y to have unit variance
y = y/ std(y);

% Generate noise signal eta
sigma = 0.1;
eta = sigma * randn(N, 1);

% Create desired signal z by adding noise to y
z = y + eta;

% Initialize arrays for storing filter coefficients
Nw = 4;      % Filter order
w_all = zeros(Nw+1, N+1);  % Matrix containing the evolution of the adaptive weights in time
e_all = zeros(N, 1);       % Error signal
% Apply LMS algorithm for mu = 0.01 and store the filter coefficients over time
mu = 0.01;
[yhat, e, w] = lms(x, z, mu, Nw);
w_all(:,1:N+1) = w;
e_all = e_all + e.^2;
% Plot results
figure;
subplot(2,1,1);
plot(0:N, w_all);
title('Evolution of filter coefficients for \mu = 0.1');
xlabel('Time');
ylabel('Filter coefficients');
legend('w[0]', 'w[1]', 'w[2]', 'w[3]', 'w[4]');
subplot(2,1,2);
plot(0:N-1, e_all);
title('Least Mean Squared Estimate Error');
xlabel('Time');
ylabel('Estimated Error');





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