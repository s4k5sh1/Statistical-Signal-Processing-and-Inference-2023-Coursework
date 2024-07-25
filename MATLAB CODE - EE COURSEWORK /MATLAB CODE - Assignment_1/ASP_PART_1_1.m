clear all
x = rand(1000,1);
figure
stem(1:1000,x)
xlabel('Sample Size (n)')
ylabel('Probability')

% expected value (theoretical mean)
theoretical_mean = (0 + 1)/2

% sample mean
x = rand(1000,1);
sample_mean = mean(x)


% theoretical standard deviation
theoretical_std = sqrt((1-0)^2/12)
% sample standard deviation
x = rand(1000,1);
sample_std = std(x)

%%

% generate ten 1000-sample realizations of X
clear all
n = 1000;
realizations = 10;
X = rand(n,realizations)

% calculate sample means  for each realization
sample_means = mean(X,1);

% calculate bias for mean
bias_mean = 0.5 - sample_means;

% plot mean error for ten realizations
figure
stem(1:realizations, bias_mean, '-o')
xlabel("Realizations of X")
ylabel("Mean Error")


% calculate sample sd for each realization

sample_sd = std(X,0,1);

% calculate bias for standard deviation

bias_sd = 0.2887 - sample_sd;

% plot the sd error for ten realizations
figure
stem(1:realizations, bias_sd, '-o')
xlabel("Realizations of X")
ylabel("Standard Deviation Error")

%%


%plot the histogram of a uniform random variable
clear all
x = rand(1000,1);
figure
histogram(x,'Normalization','probability')%normalised by number of samples taken
hold on 
a = 0; % lower bound of uniform distribution
b = 1; % upper bound of uniform distribution
x = linspace(a, b, 1000); % points at which to evaluate the pdf
% calculate the pdf
pdf = (ones(1, 1000)/(b-a))/10;
plot(x, pdf, 'LineWidth', 2, 'Color', 'r');
xlabel("x")
ylabel("Probability Density")


%%
%analysis using gaussian random variable 
clear all
x = randn(1000,1);
figure
stem(1:1000,x)
xlabel('Sample Size (n)')
ylabel('Probability')


% sample mean
x = randn(1000,1);
sample_mean = mean(x)

% sample standard deviation
x = randn(1000,1);
sample_std = std(x)

% generate ten 1000-sample realizations of X
n = 1000;
realizations = 10;
X = randn(n,realizations)

% calculate sample means  for each realization
sample_means = mean(X,1);

% calculate bias for mean
bias_mean = 0 - sample_means;

% calculate sample sd for each realization

sample_sd = std(X,0,1);

% calculate bias for standard deviation

bias_sd = 1 - sample_sd;



% plot mean error for ten realizations
figure
subplot(2,1,1)
stem(1:realizations, bias_mean, '-o')
xlabel("Realizations of X")
ylabel("Mean Error")
hold on 
subplot(2,1,2)
% plot the sd error for ten realizations
stem(1:realizations, bias_sd, '-o')
xlabel("Realizations of X")
ylabel("Standard Deviation Error")
hold off

%%

% Generate 10 realizations of 1000 samples of Gaussian random variable
clear all
mu = 0;
sigma = 1;
X = randn(1000,1);
% Plot the histogram of the samples
figure
histogram(X, 'Normalization', 'probability')

% Add the theoretical Gaussian PDF
hold on
x = linspace(mu-3*sigma, mu+3*sigma, 100);
pdf = (normpdf(x, mu, sigma));
plot(x, pdf, 'r', 'LineWidth', 2)
xlabel("x")
ylabel("Probability Desnity Function")






















