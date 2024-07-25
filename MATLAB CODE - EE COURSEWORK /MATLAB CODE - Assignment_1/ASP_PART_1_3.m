%% estimating pdf - stationary 
% function called pdf is created 
X = randn(1,100000); % Generate N (100000) samples of a Gaussian random variable
[counts, bins] = pdf(X); % pdf function to estimate the pdf of the data

%% estimating pdf - stationary and ergodic

clear all
M = 1;
N = 100;
v = rp3(M, N);
subplot(4,1,1)
[counts, bins] = pdf(v);
xlabel("x")
ylabel('Probability Density')
title("N = 100")
M = 1;
N = 1000;
v = rp3(M, N);
subplot(4,1,2)
[counts, bins] = pdf(v);
xlabel("x")
ylabel('Probability Density')
title("N = 1000")
M = 1;
N = 10000;
v = rp3(M, N);
subplot(4,1,3)
[counts, bins] = pdf(v);
xlabel("x")
ylabel('Probability Density')
title("N = 10000")

a = 0.5;
m = 3;
% Define the range of values for the x-axis
x = linspace(-0.5, 3.5, 1000);

% Calculate the uniform pdf for the defined range
y = ones(1, 1000)/m;

% Plot the theoretical pdf of the uniform random variable
subplot(4,1,4)
plot(x, y);
xlabel('x');
ylabel('Probability Density');
title('Theoretical PDF of RP3');

%% estimating pdf of non-stationary
clear all
M = 1;
N = 5000;
v = rp2(M, N);
figure;
[counts, bins] = pdf(v);
xlabel("x")
ylabel('Probability Density')
title("N = 10000")


%% functions

 function v=rp1(M,N);
a=0.02;
b=5;
Mc=ones(M,1)*b*sin((1:N)*pi/N);
Ac=a*ones(M,1)*[1:N];
v=(rand(M,N)-0.5).*Mc+Ac; 

 end


function v=rp2(M,N)
Ar=rand(M,1)*ones(1,N);
Mr=rand(M,1)*ones(1,N);
v=(rand(M,N)-0.5).*Mr+Ar;
 
end

 function v=rp3(M,N)
a=0.5;
m=3;
v=(rand(M,N)-0.5)*m + a
 end 

 function [counts, bins] = pdf(v)
    N = length(v);
    [counts, bins] = hist(v, 100); % Compute histogram of data with 100 bins
    counts = counts / trapz(bins,counts); % Normalize the histogram by number of samples considered 
    bar(bins, counts); % Plot histogram bars
    xlabel('x');
    ylabel('Probability');
end