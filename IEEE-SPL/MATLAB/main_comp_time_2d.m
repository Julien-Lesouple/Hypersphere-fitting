clear
close all
clc

rng(2)

%% Parameters
% Dimension
n = 100; % observations
d = 2; % features

%Noise
sigma = 2;

Nmc=500;

n_iter = 40;

time_em = zeros(1,Nmc);
time_kasa = zeros(1,Nmc);
time_elandau = zeros(1,Nmc);
time_iml = zeros(1,Nmc);

for ind_mc = 1:Nmc
    % Ground truth sphere
    c= randi(10,1,d)-5;
    r = randi([1 10],1,1);%

    %% Generate data
    % Hidden parameters
    kappa=0;
    mu=pi/4;
    mu_vec = [cos(mu), sin(mu)];
    u=vmrand(mu, kappa, n, 1);
    pt=[cos(u) sin(u)];

    z = c + r*pt;

    % Observations
    a = z+sigma*randn(n,d);

    %% EM
    tic
    thetas_em = EM_VmF(a,kappa,mu_vec, n_iter);
    time_em(ind_mc) = toc;

    %% E-Landau
    tic
    thetas_elandau = e_landau(a);
    time_elandau(ind_mc) = toc;
    
    %% Kasa
    tic
    thetas_kasa = kasa(a);
    time_kasa(ind_mc) = toc;

    %% IML
    tic
    thetas_iml = iml(a, n_iter);
    time_iml(ind_mc) = toc;
end

mean_time_em = mean(time_em);
mean_time_elandau = mean(time_elandau);
mean_time_kasa = mean(time_kasa);
mean_time_iml = mean(time_iml);