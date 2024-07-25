% 4.4.1 and % 4.4.2 


N = 1000; % Number of samples
a = [1, 0.9, 0.2]; % AR model coefficients
mu = 0.01; % LMS step size
Nw = 2; % Filter order (2, as we want to estimate a1 and a2)

% Generate white noise
wn = randn(N, 1);


signal = filter(1, a, wn);
lms_modes = {'signed error', 'signed regressor', 'sign sign'};

% Apply the LMS algorithm
[yhat, e, w] = lms(signal(1:end-1), signal(2:end), mu, Nw);



% Plot the convergence of the estimated coefficients
figure;
plot(w(1, :));
hold on;
plot(w(2, :));
title('Convergence of Estimated Coefficients for \mu = 0.01');
xlabel('Time');
ylabel('Coefficient Value');
legend('a1', 'a2');






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
