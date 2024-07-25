% ACF of unbiased, 1000 sample white gaussian noise

x = randn(1000,1);
[r,lags]= xcorr(x,'unbiased');
figure(1);
ACF_plot = plot(lags,r);
set(ACF_plot, 'Marker', 'none')
title('Unbiased estimate of ACF for 1000 WGN samples')
ylabel('Correlation')
xlabel('τ (Correlation lag)')
axis([-999 999 -0.8 1.2])

%%
% using zoom function to focus onto the region |τ | < 50

figure(2);
ACF_plot = plot(lags,r);
set(ACF_plot, 'Marker', 'none')
title('Unbiased estimate of ACF for 1000 WGN samples')
ylabel('Correlation')
xlabel('τ (Correlation lag)')
% Zoom in on the x-axis
zoom on;
axis([-50, 50, ylim]);
zoom off;
%%
% Generate a 1000-sample WGN vector x and filter it by a moving average (MA) filter with unit coefficients of [10]
% of order 9 and plot the acf

clear all;
x = randn(1000,1);
y=filter(ones(9,1),[1],x);
[r,lags]= xcorr(y,'unbiased');
figure(3);
ACF = stem(lags,r);
title('ACF estimate for 9th order MA filter applied on 1000 samples of WGN')
ylabel('Correlation')
xlabel('τ (Correlation lag)')
xlim([-20 20])