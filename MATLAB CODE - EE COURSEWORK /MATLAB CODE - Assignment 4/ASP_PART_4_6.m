%4.6
% Parameters
N = 1000; % Number of samples
a = [1, 0.9, 0.2]; % AR model coefficients
mu = 0.03; % LMS step size
Nw = 2; % Filter order (2, as we want to estimate a1 and a2)

% Generate white noise
wn = randn(N, 1);

% Synthesize the signal using the AR model
signal = filter(1, a, wn);
lms_modes = {'signed error', 'signed regressor', 'sign sign'};

% Apply the LMS algorithm
[yhat, e, w] = sign_lms(signal(1:end-1), signal(2:end), mu, Nw,lms_modes{3});

% Display the final estimated coefficients
final_w = w(:, end);

% Plot the convergence of the estimated coefficients
figure;
plot(w(1, :));
hold on;
plot(w(2, :));
title('Convergence of Estimated Coefficients for \mu = 0.01 using Sign-Sign');
xlabel('Time');
ylabel('Coefficient Value');
legend('a1', 'a2');


% 4.6 on audio signals 

% Read the audio file for sound "e"
[x, Fs] = audioread('sound_e.wav');
x = x(1:N);
x = x(1:N) / max(abs(x(1:N)));

% Generate desired signal z(n)
z = x(2:end);

% Adaptation gain
mu = 0.02;

% Basic LMS algorithm
[yhat_lms, e_lms, w_lms] = lms(x(1:end-1), z, mu, 2);

% Sign LMS algorithms
[yhat_signed_error, e_signed_error, w_signed_error] = sign_lms(x(1:end-1), z, mu, 2, 'signed error');
[yhat_signed_regressor, e_signed_regressor, w_signed_regressor] = sign_lms(x(1:end-1), z, mu, 2, 'signed regressor');
[yhat_sign_sign, e_sign_sign, w_sign_sign] = sign_lms(x(1:end-1), z, mu, 2, 'sign sign');

% Plot adaptive weights
figure;
subplot(4, 1, 1);
plot(w_lms(2:end, :)');
title('Basic LMS');
xlabel('Time ');
ylabel('Coefficients');

subplot(4, 1, 2);
plot(w_signed_error(2:end, :)');
title('Signed Error LMS');
xlabel('Time ');
ylabel('Coefficients');

subplot(4, 1, 3);
plot(w_signed_regressor(2:end, :)');
title('Signed Regressor LMS');
xlabel('Time ');
ylabel('Coefficients');

subplot(4, 1, 4);
plot(w_sign_sign(2:end, :)');
title('Sign-Sign LMS');
xlabel('Time ');
ylabel('Coefficients');





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






%% sign lms function

function [yhat, e, w] = sign_lms(x, z, mu, Nw, mode)
    % Initialization
    N = length(x);          % Length of input signal
    w = zeros(Nw+1, N);     % Matrix containing the evolution of the adaptive weights in time
    yhat = zeros(N, 1);     % Estimated output signal
    e = zeros(N, 1);        % Error signal

    % Sign LMS algorithm
    for n = 1:N
        xvec = [flipud(x(max(n-Nw,1):n)); zeros(Nw-n+1, 1)]; % Construct input vector
        yhat(n) = w(:,n)' * xvec;  % Compute estimated output
        e(n) = z(n) - yhat(n);     % Compute error
        
        % Update filter coefficients based on the selected mode
        switch mode
            case 'signed error'
                w(:,n+1) = w(:,n) + mu * sign(e(n)) * xvec;
            case 'signed regressor'
                w(:,n+1) = w(:,n) + mu * e(n) * sign(xvec);
            case 'sign sign'
                w(:,n+1) = w(:,n) + mu * sign(e(n)) * sign(xvec);
            otherwise
                error('Invalid mode. Use "signed_error", "signed_regressor", or "sign_sign".');
        end
    end
end

%%
