%% figure 1: auc vs. alpha+beta+gamma (3D) and alpha-beta/alpha-gamma/beta-gamma (2D)
hfig = figure;
subplot(2,2,1); % alpha-beta-gamma
diameter = 100*auc_f(:,tot_query_times);
color = auc_f(:,tot_query_times);
scatter3(log_alpha,beta,log_gamma,diameter,color,'linewidth',3);
xlabel('log(\alpha)');
ylabel('\beta%');
zlabel('log(\gamma)');
color_bottom = min(color);
color_top = max(color);
if color_bottom==color_top
    color_bottom = 0;
    color_top = 1;
end
caxis manual
caxis([color_bottom color_top]);
colorbar;
title(sprintf('%s',feedback_method_name));

log_alpha_set = unique(log_alpha);
beta_set = unique(beta);
log_gamma_set = unique(log_gamma);

subplot(2,2,2); % alpha-beta
alpha_beta_auc = zeros(length(log_alpha_set), length(beta_set));
for i=1:length(log_alpha_set)
    for j=1:length(beta_set)
        alpha_beta_auc(i,j) = mean(auc_f(log_alpha==log_alpha_set(i) & beta==beta_set(j),1));
    end
end
h = bar3(alpha_beta_auc');
for i = 1:length(h)
     zdata = get(h(i),'Zdata');
     set(h(i),'Cdata',zdata)
end
caxis manual
caxis([color_bottom color_top]);
colorbar;

xlabel('log(\alpha)');ylabel('\beta%');zlabel('auc');
set(gca,'xticklabel',num2cell(log_alpha_set));
set(gca,'yticklabel',num2cell(beta_set));
title('\alpha-\beta%-auc'); 

subplot(2,2,3); % alpha-gamma
alpha_gamma_auc = zeros(length(log_alpha_set), length(log_gamma_set));
for i=1:length(log_alpha_set)
    for j=1:length(log_gamma_set)
        alpha_gamma_auc(i,j) = mean(auc_f(log_alpha==log_alpha_set(i) & log_gamma==log_gamma_set(j),1));
    end
end
h = bar3(alpha_gamma_auc');
for i = 1:length(h)
     zdata = get(h(i),'Zdata');
     set(h(i),'Cdata',zdata)
end
caxis manual
caxis([color_bottom color_top]);
colorbar;
xlabel('log(\alpha)');ylabel('log(\gamma)');zlabel('auc');
set(gca,'xticklabel',num2cell(log_alpha_set));
set(gca,'yticklabel',num2cell(log_gamma_set));
title('\alpha-\gamma-auc');

subplot(2,2,4); % beta-gamma
beta_gamma_auc = zeros(length(beta_set), length(log_gamma_set));
for i=1:length(beta_set)
    for j=1:length(log_gamma_set)
        beta_gamma_auc(i,j) = mean(auc_f(beta==beta_set(i) & log_gamma==log_gamma_set(j),1));
    end
end
h = bar3(beta_gamma_auc');
for i = 1:length(h)
     zdata = get(h(i),'Zdata');
     set(h(i),'Cdata',zdata)
end
caxis manual
caxis([color_bottom color_top]);
colorbar;
xlabel('\beta%');ylabel('log(\gamma)');zlabel('auc');
set(gca,'xticklabel',num2cell(beta_set));
set(gca,'yticklabel',num2cell(log_gamma_set));
title('\beta%-\gamma-auc');

saveas(hfig, [result_mat_file(1:end-4), '-3D.fig']);
if ~show_figure_flag, close(hfig); end