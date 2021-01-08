clear
close all
clc

rng(2)

%% Parameters
% Dimensions
n = 100; % observations
d = 2; % features

% Prior
kappa=2; % 0: uniform prior, >0: Von Mises Fisher, in the article: kappa = 2
mu=pi/4;
mu_vec = [cos(mu), sin(mu)];

%Noise
sigma2_vec = 10.^(-2:0.1:1);

% Monte-Carlo
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

err_r_elandau = zeros(1,length(sigma2_vec));
err_c_elandau = zeros(1,length(sigma2_vec));
err_rc_elandau = zeros(1,length(sigma2_vec));

std_r_elandau = zeros(1,length(sigma2_vec));
std_c_elandau = zeros(1,length(sigma2_vec));
std_rc_elandau = zeros(1,length(sigma2_vec));

err_r_kasa = zeros(1,length(sigma2_vec));
err_c_kasa = zeros(1,length(sigma2_vec));
err_rc_kasa = zeros(1,length(sigma2_vec));

std_r_kasa = zeros(1,length(sigma2_vec));
std_c_kasa = zeros(1,length(sigma2_vec));
std_rc_kasa = zeros(1,length(sigma2_vec));

err_r_iml = zeros(1,length(sigma2_vec));
err_c_iml = zeros(1,length(sigma2_vec));
err_rc_iml = zeros(1,length(sigma2_vec));

std_r_iml = zeros(1,length(sigma2_vec));
std_c_iml = zeros(1,length(sigma2_vec));
std_rc_iml = zeros(1,length(sigma2_vec));

%% Loop over sigma
for ind_sigma = 1:length(sigma2_vec)
    sigma = sqrt(sigma2_vec(ind_sigma));
    disp(['sigma value ', num2str(ind_sigma),' over ',num2str(length(sigma2_vec))])
    thetas_diff_em = zeros(d+2,Nmc);
    thetas_diff_kasa = zeros(d+1,Nmc);
    thetas_diff_elandau = zeros(d+1,Nmc);
    thetas_diff_iml = zeros(d+1,Nmc);
    
    % Monte-Carlo
    for ind_mc = 1:Nmc
        %% Generate data
        % Ground truth sphere
        c= randi(10,1,d)-5;
        r = randi([1 10],1,1);
        
        % Hidden parameters
        u=vmrand(mu, kappa, n, 1);
        pt=[cos(u) sin(u)];

        z = c + r*pt;

        % Observations
        a = z+sigma*randn(n,d);
        
        % Parameters to estimate
        theta_th = [r;c';sigma^2];

        %% EM
        theta_em = EM_VmF(a,kappa,mu_vec, n_iter);
        thetas_diff_em(:,ind_mc) = theta_th-theta_em(:,end);
        %% E-Landau
        thetas_diff_elandau(:,ind_mc) = theta_th(1:end-1)-e_landau(a);

        %% Kasa
        thetas_diff_kasa(:,ind_mc) = theta_th(1:end-1)-kasa(a);

        %% IML
        theta_iml = iml(a, n_iter);
        thetas_diff_iml(:,ind_mc) = theta_th(1:end-1)-theta_iml(:,end);
    end
    
    %% Compute Monte-Carlo mean and standard deviations
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

    % E-Landau
    err_r_elandau(ind_sigma) = mean(vecnorm(thetas_diff_elandau(1,:),2,1).^2);
    err_c_elandau(ind_sigma) = mean(vecnorm(thetas_diff_elandau(2:d+1,:),2,1).^2);
    err_rc_elandau(ind_sigma) = mean(vecnorm(thetas_diff_elandau(1:d+1,:),2,1).^2);

    std_r_elandau(ind_sigma) = std(vecnorm(thetas_diff_elandau(1,:),2,1).^2);
    std_c_elandau(ind_sigma) = std(vecnorm(thetas_diff_elandau(2:d+1,:),2,1).^2);
    std_rc_elandau(ind_sigma) = std(vecnorm(thetas_diff_elandau(1:d+1,:),2,1).^2);

    % Kasa
    err_r_kasa(ind_sigma) = mean(vecnorm(thetas_diff_kasa(1,:),2,1).^2);
    err_c_kasa(ind_sigma) = mean(vecnorm(thetas_diff_kasa(2:d+1,:),2,1).^2);
    err_rc_kasa(ind_sigma) = mean(vecnorm(thetas_diff_kasa(1:d+1,:),2,1).^2);

    std_r_kasa(ind_sigma) = std(vecnorm(thetas_diff_kasa(1,:),2,1).^2);
    std_c_kasa(ind_sigma) = std(vecnorm(thetas_diff_kasa(2:d+1,:),2,1).^2);
    std_rc_kasa(ind_sigma) = std(vecnorm(thetas_diff_kasa(1:d+1,:),2,1).^2);

    % IML
    err_r_iml(ind_sigma) = mean(vecnorm(thetas_diff_iml(1,:),2,1).^2);
    err_c_iml(ind_sigma) = mean(vecnorm(thetas_diff_iml(2:d+1,:),2,1).^2);
    err_rc_iml(ind_sigma) = mean(vecnorm(thetas_diff_iml(1:d+1,:),2,1).^2);

    std_r_iml(ind_sigma) = std(vecnorm(thetas_diff_iml(1,:),2,1).^2);
    std_c_iml(ind_sigma) = std(vecnorm(thetas_diff_iml(2:d+1,:),2,1).^2);
    std_rc_iml(ind_sigma) = std(vecnorm(thetas_diff_iml(1:d+1,:),2,1).^2);
end

% plot
plot_figures_2D
