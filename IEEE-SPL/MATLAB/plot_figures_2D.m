close all

plot_ftsize = 21;

%% Figures in article
figure
hold on
plot(10*log10(sigma2_vec), 10*log10(err_rc_em), 'linewidth',3)
plot(10*log10(sigma2_vec), 10*log10(err_rc_iml), '--', 'linewidth',3)
plot(10*log10(sigma2_vec), 10*log10(err_rc_elandau), ':', 'linewidth',3)
plot(10*log10(sigma2_vec), 10*log10(err_rc_kasa),'-.', 'linewidth',3)

axis on
box on
grid on
axis square

set(gca, 'fontsize', plot_ftsize, 'linewidth', 1.5,'GridLineStyle','--')
xlabel('$10\log_{10}{(\sigma^2)}$','FontSize',plot_ftsize,'interpreter','latex')
ylabel('$10\log_{10}{\mbox{MSE}}$','FontSize',plot_ftsize,'interpreter','latex')
title('MSE for ($r$,{\boldmath$c$}) in 2D','FontSize',plot_ftsize,'interpreter','latex')
legend({'EM', 'IML', 'E-Landau', 'Kasa'},'Fontsize',plot_ftsize,'interpreter','latex')


%% Figures in technical report

% Radius
figure
hold on
plot(10*log10(sigma2_vec), 10*log10(err_r_em), 'linewidth',3)
plot(10*log10(sigma2_vec), 10*log10(err_r_iml), '--', 'linewidth',3)
plot(10*log10(sigma2_vec), 10*log10(err_r_elandau), ':', 'linewidth',3)
plot(10*log10(sigma2_vec), 10*log10(err_r_kasa),'-.', 'linewidth',3)

axis on
box on
grid on
axis square

set(gca, 'fontsize', plot_ftsize, 'linewidth', 1.5,'GridLineStyle','--')
xlabel('$10\log_{10}{(\sigma^2)}$','FontSize',plot_ftsize,'interpreter','latex')
ylabel('$10\log_{10}{\mbox{MSE}}$','FontSize',plot_ftsize,'interpreter','latex')
title('MSE for $r$ in 2D','FontSize',plot_ftsize,'interpreter','latex')
legend({'EM', 'IML', 'E-Landau', 'Kasa'},'Fontsize',plot_ftsize,'interpreter','latex')

% Center
figure
hold on
plot(10*log10(sigma2_vec), 10*log10(err_c_em), 'linewidth',3)
plot(10*log10(sigma2_vec), 10*log10(err_c_iml), '--', 'linewidth',3)
plot(10*log10(sigma2_vec), 10*log10(err_c_elandau), ':', 'linewidth',3)
plot(10*log10(sigma2_vec), 10*log10(err_c_kasa),'-.', 'linewidth',3)

axis on
box on
grid on
axis square

set(gca, 'fontsize', plot_ftsize, 'linewidth', 1.5,'GridLineStyle','--')
xlabel('$10\log_{10}{(\sigma^2)}$','FontSize',plot_ftsize,'interpreter','latex')
ylabel('$10\log_{10}{\mbox{MSE}}$','FontSize',plot_ftsize,'interpreter','latex')
title('MSE for {\boldmath$c$} in 2D','FontSize',plot_ftsize,'interpreter','latex')
legend({'EM', 'IML', 'E-Landau', 'Kasa'},'Fontsize',plot_ftsize,'interpreter','latex')


% Variance
figure
hold on
plot(10*log10(sigma2_vec), 10*log10(err_sigma2_em), 'linewidth',3)

axis on
box on
grid on
axis square
set(gca, 'fontsize', plot_ftsize, 'linewidth', 1.5,'GridLineStyle','--')
xlabel('$10\log_{10}{(\sigma^2)}$','FontSize',plot_ftsize,'interpreter','latex')
ylabel('$10\log_{10}{\mbox{MSE}}$','FontSize',plot_ftsize,'interpreter','latex')
title('MSE for $\sigma^2$ in 2D','FontSize',plot_ftsize,'interpreter','latex')
