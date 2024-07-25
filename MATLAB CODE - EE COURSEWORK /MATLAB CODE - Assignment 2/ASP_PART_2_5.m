% 2.5
%loading all csv files 

clear all;
load ECG1_Matteo.csv
load ECG2_Matteo.csv
load ECG3_Matteo.csv

%%
% plotting the whole ECG data
ECG_data= ECG1_Matteo(:, 3);


%%
% splitting the data
unconstrained_breathing = ECG_data(8137:126359);
constrained_50 = ECG_data(130111:251349);
constrained_15 = ECG_data(254857:376555);

fsECG = 500
[xRRI,fsRRI]=ECG_to_RRI(constrained_15,fsECG);


%%
rr_unconstrained = xRRI;
% rr_constrained_50 = xRRI_2;
% rr_constrained_15 = xRRI_3;

h = 60 ./rr_unconstrained ; % calculates heart rate 
figure;
histogram(h, 'Normalization','pdf');
title('PDE of Original Heart Rate');
ylabel('Probability');
xlabel('Samples');


%%

% Calculate the averaged heart rates with different alpha values
alpha1 = 1;
alpha2 = 0.6;
n = 10;

h_avg1 = average_hr(h, n, alpha1);
h_avg2 = average_hr(h, n, alpha2);


figure;
subplot(3,1,1);
histogram(h, 'Normalization', 'pdf');
title('PDE of Original Heart Rate');
ylabel('Probability');
xlabel('Samples');

subplot(3,1,2);
histogram(h_avg1, 'Normalization', 'pdf','FaceColor',"#D95319");
title('PDE of Averaged Heart Rate with \alpha = 1');
ylabel('Probability');
xlabel('Samples');

subplot(3,1,3);
histogram(h_avg2, 'Normalization', 'pdf','FaceColor',"#D95319");
title('PDE of Averaged Heart Rate with \alpha = 0.6');
ylabel('Probability');
xlabel('Samples');






%%
% autocorrelation of RRI data 
rr_unconstrained_detrend = detrend(rr_unconstrained);
% rr_constrained_50_detrend = detrend(rr_constrained_50);
% rr_constrained_15_detrend = detrend(rr_constrained_15);
[autocorr_rr_unconstrained_detrend , lags1] = xcorr(rr_unconstrained_detrend, 'coeff');
figure;
plot(lags1, autocorr_rr_unconstrained_detrend);
title('Autocorrelation of Constrained Breathing Data ');
xlabel('Lag');
ylabel('Autocorrelation');


%%
% finding ar model of heart rate trials 

max_p = 10; 
mdl = zeros(max_p, 1);
aic = zeros(max_p, 1);
aicc = zeros(max_p, 1);

for p = 1:max_p
    
    [a, e, k] = aryule(rr_unconstrained_detrend, p);
    
    
    n = length(rr_unconstrained_detrend);
    mdl(p) = log(e) + p*log(n)/n;
    aic(p) = log(e) + 2*p/n;
    aicc(p) = aic(p) + (2*p*(p+1))/(n-p-1);
end


figure;
subplot(2,1,1)
plot(1:max_p, mdl, 'o-', 1:max_p, aic, 'o-', 1:max_p, aicc, 'o-');
legend('MDL', 'AIC', 'AICc');
title('MDL,AIC and AICc from p = 0 to p = 10')
xlabel('Model Order');
ylabel('Critical Value');


x= rr_unconstrained;
[arcoefs,E,K] = aryule(x,10);
pacf = -K;
subplot(2,1,2)
stem(pacf)
xlabel('Lag')
ylabel('Partial ACF')
title('Partial Autocorrelation Sequence')
xlim([1 10])
hold on
x_std = rr_unconstrained_detrend;
[arcoefs,E,K] = aryule(x_std,10);
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
function avg_h = average_hr(h, n, alpha)
    avg_h = zeros(length(h) - n + 1, 1);
    for i = 1:(length(h) - n + 1)
        avg_h(i) = (1/n) * alpha * sum(h(i:(i+n-1)));
    end
end


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