% Estimate the CCF for the sequences x and y generated above 

clear all;
x = randn(1000,1);
y=filter(ones(9,1),[1],x);
[r,lags]= xcorr(x,y,'unbiased');
figure(4);
ACF = stem(lags,r);
title('CCF estimate for 9th order MA filter applied on 1000 samples of WGN')
ylabel('Correlation')
xlabel('τ (Correlation lag)')
xlim([-20 20])