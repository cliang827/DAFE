%% figure 1: auc vs. alpha+beta+gamma (3D) and alpha-beta/alpha-gamma/beta-gamma (2D)
hfig = figure;
subplot(2,2,1); % alpha-beta-gamma
diameter = 100*auc_f_kmean(:,tot_query_times);
color = auc_f_kmean(:,tot_query_times);
scatter3(alpha,beta,gamma,diameter,color,'linewidth',3);
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
title(sprintf('%s+%d',feedback_method,v_sum_constraint));

alpha_set = unique(alpha);
beta_set = unique(beta);
gamma_set = unique(gamma);

subplot(2,2,2); % alpha-beta
alpha_beta_auc = zeros(length(alpha_set), length(beta_set));
for i=1:length(alpha_set)
    for j=1:length(beta_set)
        alpha_beta_auc(i,j) = mean(auc_f_kmean(alpha==alpha_set(i) & beta==beta_set(j),2));
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
set(gca,'xticklabel',num2cell(alpha_set));
set(gca,'yticklabel',num2cell(beta_set));
title('\alpha-\beta%-auc'); 

subplot(2,2,3); % alpha-gamma
alpha_gamma_auc = zeros(length(alpha_set), length(gamma_set));
for i=1:length(alpha_set)
    for j=1:length(gamma_set)
        alpha_gamma_auc(i,j) = mean(auc_f_kmean(alpha==alpha_set(i) & gamma==gamma_set(j),2));
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
set(gca,'xticklabel',num2cell(alpha_set));
set(gca,'yticklabel',num2cell(gamma_set));
title('\alpha-\gamma-auc');

subplot(2,2,4); % beta-gamma
beta_gamma_auc = zeros(length(beta_set), length(gamma_set));
for i=1:length(beta_set)
    for j=1:length(gamma_set)
        beta_gamma_auc(i,j) = mean(auc_f_kmean(beta==beta_set(i) & gamma==gamma_set(j),2));
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
set(gca,'yticklabel',num2cell(gamma_set));
title('\beta%-\gamma-auc');

saveas(hfig, [result_mat_file(1:end-4), '-3D.fig']);
if ~show_figure_flag, close(hfig); end