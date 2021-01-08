clear
close all
clc

rng(2)

%% Parameters
% Dimension
n = 100; % observations
d = 3; % features

%Noise
sigma = 2;

Nmc=500;

n_iter = 40;

time_em = zeros(1,Nmc);
time_ils = zeros(1,Nmc);
time_fgfa = zeros(1,Nmc);

for ind_mc = 1:Nmc
    % Ground truth sphere
    c= randi(10,1,d)-5;
    r = randi([1 10],1,1);%

    %% Generate data
    % Hidden parameters
    kappa=0;
    mu=pi/4;
    mu_vec = randn(1,3);
    mu_vec = mu_vec./norm(mu_vec);
    u=randn(n,d);
    pt=u./vecnorm(u,2,2);

    z = c + r*pt;

    % Observations
    a = z+sigma*randn(n,d);

    %% EM
    tic
    thetas_em = EM_VmF(a,kappa,mu_vec, n_iter);
    time_em(ind_mc) = toc;

    %% ILS
    tic
    thetas_ils = iml(a, n_iter);
    time_ils(ind_mc) = toc;
    
    %% FGFA
    tic
    thetas_fgfa = fgfa(a);
    time_fgfa(ind_mc) = toc;
end

mean_time_em = mean(time_em);
mean_time_ils = mean(time_ils);
mean_time_fgfa = mean(time_fgfa);
