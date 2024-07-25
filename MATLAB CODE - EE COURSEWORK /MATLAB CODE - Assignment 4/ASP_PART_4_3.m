% 4.3 - gear shifting 

sd = 0.1;
adap_gain = 0.01;

a = 1;
b = [1 2 3 2 1];

N = 1000;
N_w = length(b);

x = randn(N,1);
y = filter(b, a, x);
y = y./sqrt(sum(b.*b));

n = sd*randn(N,1);
z = y+n;

[y_e,e,w] = lms_gs(x, z, adap_gain, N_w+1);
w_sc = w.*sqrt(sum(b.*b));
error_squared = e.^2;

figure;
hold on
plot(w_sc(1,:), 'DisplayName', 'w[0]')
plot(w_sc(2,:), 'DisplayName', 'w[1]')
plot(w_sc(3,:), 'DisplayName', 'w[2]')
plot(w_sc(4,:), 'DisplayName', 'w[3]')
plot(w_sc(5,:), 'DisplayName', 'w[4]')
title('Time evolution of coefficients')
xlabel('Sample')
ylabel('Coefficients')
legend('show')
hold off;