% 2.4.1 
clear all;
load NASDAQ.mat; 
% Standardize the data to zero mean and unit variance
NASDAQ_std = (NASDAQ.Close - mean(NASDAQ.Close)) / std(NASDAQ.Close);


max_p = 10; % Maximum model order
mdl = zeros(max_p, 1);
aic = zeros(max_p, 1);
aicc = zeros(max_p, 1);

for p = 1:max_p
 
    [a, e, k] = aryule(NASDAQ_std, p);
    

    n = length(NASDAQ_std);
    mdl(p) = log(e) + p*log(n)/n;
    aic(p) = log(e) + 2*p/n;
    aicc(p) = aic(p) + (2*p*(p+1))/(n-p-1);
end

figure;
subplot(1,2,1)
plot(1:max_p, mdl, 'o-', 1:max_p, aic, 'o-', 1:max_p, aicc, 'o-');
legend('MDL', 'AIC', 'AICc');
title('MDL,AIC and AICc from p = 0 to p = 10')
xlabel('Model Order');
ylabel('Critical Value');


x = NASDAQ.Close; % Extract data from second column
[arcoefs,E,K] = aryule(x,15);
pacf = -K;
subplot(1,2,2)
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

%%

%%