
%% rp1 - stationarity 
clear all
M = 100;
N = 100;
v = rp1(M, N);
mean_rp1 = mean(v,1); % Ensemble mean
std_rp1 = std(v,0,1); % Ensemble standard deviation
time = 1:N; % Time

% Plotting
figure
subplot(1,2,1)
plot(time, mean_rp1,'Color', 'b')
xlabel('Time')
ylabel('Ensemble mean')
subplot(1,2,2)
plot(time, std_rp1, 'Color', 'b')
xlabel('Time')
ylabel('Ensemble standard deviation')




 %% rp2 - stationarity 

clear all

M = 100;
N = 100;
v = rp2(M, N);
mean_rp2 = mean(v,1); % Ensemble mean
std_rp2 = std(v,0,1); % Ensemble standard deviation
time = 1:N; % Time

% Plotting
figure
subplot(1,2,1)
plot(time, mean_rp2,'Color', 'b')
xlabel('Time')
ylabel('Ensemble mean')
subplot(1,2,2)
plot(time, std_rp2, 'Color', 'b')
xlabel('Time')
ylabel('Ensemble standard deviation')

average_over_time = mean(mean_rp2)
sd_over_time = mean(std_rp2)



%% rp3 - stationarity 

clear all

M = 100;
N = 100;
v = rp3(M, N);
mean_rp3 = mean(v,1); % Ensemble mean
std_rp3 = std(v,0,1); % Ensemble standard deviation
time = 1:N; % Time

% Plotting
figure
subplot(1,2,1)
plot(time, mean_rp3,'Color', 'b')
xlabel('Time')
ylabel('Ensemble mean')
subplot(1,2,2)
plot(time, std_rp3, 'Color', 'b')
xlabel('Time')
ylabel('Ensemble standard deviation')

average_over_time = mean(mean_rp3)
sd_over_time = mean(std_rp3)



%% rp1 - ergodicity 
clear all
M = 4;
N = 1000;
v = rp1(M, N);
for i = 1:M
    mean_rp1(i) = mean(v(i,:)) % Calculate mean for each realization
    std_rp1(i) = std(v(i,:)) % Calculate standard deviation for each realization
end

mean_mean_total = mean(mean_rp1); % average mean of all realizations
mean_std_total = mean(std_rp1); % average standard deviation of all realizatio


%% rp2 - ergodicity 
clear all
M = 4;
N = 1000;
v = rp2(M, N);
for i = 1:M
    mean_rp2(i) = mean(v(i,:)) % Calculate mean for each realization
    std_rp2(i) = std(v(i,:)) % Calculate standard deviation for each realization
end

mean_mean_total = mean(mean_rp2); % average mean of all realizations
mean_std_total = mean(std_rp2); % average standard deviation of all realizatio


%% rp3 - ergodicity 
clear all
M = 4;
N = 1000;
v = rp3(M, N);
for i = 1:M
    mean_rp3(i) = mean(v(i,:)) % Calculate mean for each realization
    std_rp3(i) = std(v(i,:)) % Calculate standard deviation for each realization
end

mean_mean_total = mean(mean_rp3); % average mean of all realizations
mean_std_total = mean(std_rp3); % average standard deviation of all realizatio

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

