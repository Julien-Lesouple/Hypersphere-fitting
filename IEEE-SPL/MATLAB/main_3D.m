clear
close all
clc

rng(2)

addpath('yuhuichen1015-SphericalDistributionsRand-224b007') % Functions for sampling von Mises Fisher 

%% Parameters
% Dimension
n = 100; % observations
d = 3; % features

% Prior
kappa = 0;
mu = pi/4;
nu = pi/3;
mu_vec = [sin(mu)*cos(nu), sin(mu)*sin(nu), cos(mu)];

% Noise
sigma2_vec = 10.^(-2:0.1:1);

% Monte-Carlo iterations
Nmc=500;

% Algorithms
n_iter = 40;

%% Initialization
err_r_em = zeros(1,length(sigma2_vec));
err_c_em = zeros(1,length(sigma2_vec));
err_sigma2_em = zeros(1,length(sigma2_vec));
err_rc_em = zeros(1,length(sigma2_vec));
err_rcsigma_em = zeros(1,length(sigma2_vec));

std_r_em = zeros(1,length(sigma2_vec));
std_c_em = zeros(1,length(sigma2_vec));
std_sigma2_em = zeros(1,length(sigma2_vec));
std_rc_em = zeros(1,length(sigma2_vec));
std_rcsigma_em = zeros(1,length(sigma2_vec));

err_r_fgfa = zeros(1,length(sigma2_vec));
err_c_fgfa = zeros(1,length(sigma2_vec));
err_rc_fgfa = zeros(1,length(sigma2_vec));

std_r_fgfa = zeros(1,length(sigma2_vec));
std_c_fgfa = zeros(1,length(sigma2_vec));
std_rc_fgfa = zeros(1,length(sigma2_vec));

err_r_ils = zeros(1,length(sigma2_vec));
err_c_ils = zeros(1,length(sigma2_vec));
err_rc_ils = zeros(1,length(sigma2_vec));

std_r_ls = zeros(1,length(sigma2_vec));
std_c_ils = zeros(1,length(sigma2_vec));
std_rc_ils = zeros(1,length(sigma2_vec));

%% Loop over sigma
for ind_sigma = 1:length(sigma2_vec)
    sigma = sqrt(sigma2_vec(ind_sigma));
    disp(['sigma value ', num2str(ind_sigma),' over ',num2str(length(sigma2_vec))])
    thetas_diff_em = zeros(d+2,Nmc);
    thetas_diff_fgfa = zeros(d+1,Nmc);
    thetas_diff_ils = zeros(d+1,Nmc);
    
    % Monte-Carlo
    for ind_mc = 1:Nmc
        %% Generate data
        % Ground truth sphere
        c= randi(10,1,d)-5;
        r = randi([1 10],1,1);
        
        % Hidden parameters
        if kappa == 0
            pt = randn(n,d);
            pt = pt./vecnorm(pt,2,2);
        else
            pt=randVMF(n, mu_vec, kappa);
        end
        
        z = c + r*pt;

        % Observations
        a = z+sigma*randn(n,d);
        
        % Parameters to estimate
        theta_th = [r;c';sigma^2];

        %% EM
        theta_em = EM_VmF(a,kappa,mu_vec, n_iter);
        thetas_diff_em(:,ind_mc) = theta_th-theta_em(:,end);
        
        %% FGFA
        thetas_diff_fgfa(:,ind_mc) = theta_th(1:end-1)-fgfa(a);

        %% ILS
        theta_ils = iml(a, n_iter);
        thetas_diff_ils(:,ind_mc) = theta_th(1:end-1)-theta_ils(:,end);
    end

    %% Compute Monte-Carlo means and standard deviations
    % EM
    err_r_em(ind_sigma) = mean(vecnorm(thetas_diff_em(1,:),2,1).^2);
    err_c_em(ind_sigma) = mean(vecnorm(thetas_diff_em(2:d+1,:),2,1).^2);
    err_sigma2_em(ind_sigma) = mean(vecnorm(thetas_diff_em(d+2,:),2,1).^2);
    err_rc_em(ind_sigma) = mean(vecnorm(thetas_diff_em(1:d+1,:),2,1).^2);
    err_rcsigma_em(ind_sigma) = mean(vecnorm(thetas_diff_em,2,1).^2);

    std_r_em(ind_sigma) = std(vecnorm(thetas_diff_em(1,:),2,1).^2);
    std_c_em(ind_sigma) = std(vecnorm(thetas_diff_em(2:d+1,:),2,1).^2);
    std_sigma2_em(ind_sigma) = std(vecnorm(thetas_diff_em(d+2,:),2,1).^2);
    std_rc_em(ind_sigma) = std(vecnorm(thetas_diff_em(1:d+1,:),2,1).^2);
    std_rcsigma_em(ind_sigma) = std(vecnorm(thetas_diff_em,2,1).^2);

    % FGFA
    err_r_fgfa(ind_sigma) = mean(vecnorm(thetas_diff_fgfa(1,:),2,1).^2);
    err_c_fgfa(ind_sigma) = mean(vecnorm(thetas_diff_fgfa(2:d+1,:),2,1).^2);
    err_rc_fgfa(ind_sigma) = mean(vecnorm(thetas_diff_fgfa(1:d+1,:),2,1).^2);

    std_r_fgfa(ind_sigma) = std(vecnorm(thetas_diff_fgfa(1,:),2,1).^2);
    std_c_fgfa(ind_sigma) = std(vecnorm(thetas_diff_fgfa(2:d+1,:),2,1).^2);
    std_rc_fgfa(ind_sigma) = std(vecnorm(thetas_diff_fgfa(1:d+1,:),2,1).^2);

    % ILS
    err_r_ils(ind_sigma) = mean(vecnorm(thetas_diff_ils(1,:),2,1).^2);
    err_c_ils(ind_sigma) = mean(vecnorm(thetas_diff_ils(2:d+1,:),2,1).^2);
    err_rc_ils(ind_sigma) = mean(vecnorm(thetas_diff_ils(1:d+1,:),2,1).^2);

    std_r_ls(ind_sigma) = std(vecnorm(thetas_diff_ils(1,:),2,1).^2);
    std_c_ils(ind_sigma) = std(vecnorm(thetas_diff_ils(2:d+1,:),2,1).^2);
    std_rc_ils(ind_sigma) = std(vecnorm(thetas_diff_ils(1:d+1,:),2,1).^2);
end

%plot
plot_figures_3D
