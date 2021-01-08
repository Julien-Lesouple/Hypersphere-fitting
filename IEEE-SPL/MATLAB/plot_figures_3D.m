close all

plot_ftsize = 21;

%% Figures in article
figure
hold on
plot(10*log10(sigma2_vec), 10*log10(err_rc_em), 'linewidth',3)
plot(10*log10(sigma2_vec), 10*log10(err_rc_ils), '--', 'linewidth',3)
plot(10*log10(sigma2_vec), 10*log10(err_rc_fgfa), ':', 'linewidth',3)

axis on
box on
grid on
axis square

set(gca, 'fontsize', plot_ftsize, 'linewidth', 1.5,'GridLineStyle','--')
xlabel('$10\log_{10}{(\sigma^2)}$','FontSize',plot_ftsize,'interpreter','latex')
ylabel('$10\log_{10}{\mbox{MSE}}$','FontSize',plot_ftsize,'interpreter','latex')
title('MSE on ($r$,{\boldmath$c$}) in 3D','FontSize',plot_ftsize,'interpreter','latex')
legend({'EM', 'ILS', 'FGFA'},'Fontsize',plot_ftsize,'interpreter','latex')


%% Figures in technical report

% Radius
figure
hold on
plot(10*log10(sigma2_vec), 10*log10(err_r_em), 'linewidth',3)
plot(10*log10(sigma2_vec), 10*log10(err_r_ils), '--', 'linewidth',3)
plot(10*log10(sigma2_vec), 10*log10(err_r_fgfa), ':', 'linewidth',3)

axis on
box on
grid on
axis square

set(gca, 'fontsize', plot_ftsize, 'linewidth', 1.5,'GridLineStyle','--')
xlabel('$10\log_{10}{(\sigma^2)}$','FontSize',plot_ftsize,'interpreter','latex')
ylabel('$10\log_{10}{\mbox{MSE}}$','FontSize',plot_ftsize,'interpreter','latex')
title('MSE on $r$ in 3D','FontSize',plot_ftsize,'interpreter','latex')
legend({'EM', 'ILS', 'FGFA'},'Fontsize',plot_ftsize,'interpreter','latex')

% Center
figure
hold on
plot(10*log10(sigma2_vec), 10*log10(err_c_em), 'linewidth',3)
plot(10*log10(sigma2_vec), 10*log10(err_c_ils), '--', 'linewidth',3)
plot(10*log10(sigma2_vec), 10*log10(err_c_fgfa), ':', 'linewidth',3)

axis on
box on
grid on
axis square

set(gca, 'fontsize', plot_ftsize, 'linewidth', 1.5,'GridLineStyle','--')
xlabel('$10\log_{10}{(\sigma^2)}$','FontSize',plot_ftsize,'interpreter','latex')
ylabel('$10\log_{10}{\mbox{MSE}}$','FontSize',plot_ftsize,'interpreter','latex')
title('MSE on {\boldmath$c$} in 3D','FontSize',plot_ftsize,'interpreter','latex')
legend({'EM', 'ILS', 'FGFA'},'Fontsize',plot_ftsize,'interpreter','latex')


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
title('MSE on $\sigma^2$ in 3D','FontSize',plot_ftsize,'interpreter','latex')
