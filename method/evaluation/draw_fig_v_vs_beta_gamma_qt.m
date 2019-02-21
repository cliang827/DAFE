%% figure 3: fix alpha, inspect beta and gamma (qt=1) 
alpha_star_ix = find(log_alpha==log_alpha(ix_star));
beta_num = length(unique(beta(alpha_star_ix)));
gamma_num = length(unique(log_gamma(alpha_star_ix)));
v_set = cell(beta_num*gamma_num, 1);
for i=1:length(alpha_star_ix)
    v_set{i} = squeeze(mean(v{i},3));
end

for qt=1:ctrl_para.exp_para.tot_query_times
    hfig = figure;
    for i=1:length(alpha_star_ix)
        subplot(beta_num,gamma_num,i); 

        v = v_set{i}(:,qt); %qt=1
        edges = 0:0.05:1;
        h = histogram(v,edges);
        [~, ixv]= sort(v, 'descend');
        xlabel(sprintf('log(gamma)=%.2f', log_gamma(alpha_star_ix(i))));
        ylabel(sprintf('beta%%=%.1f%%', 100*beta(alpha_star_ix(i))));

        l = legend(h, sprintf('[%.2f-%.2f:%.2f%% ]', ...
            mean(v), std(v), 100*(auc_f_kmean(alpha_star_ix(i),qt)-auc_y_kmean(alpha_star_ix(i),1))), ...
            'location', 'southeast');

        if alpha_star_ix(i)==ix_star
            title(l, 'best-parameters!', 'color', 'red');
        end
    end
    suptitle(sprintf('qt=%d with best alpha (log(alpha)=%.2f)', qt, log10(log_alpha(ix_star))));
    saveas(hfig, sprintf('%s-beta-gamma-qt%d.fig',result_mat_file(1:end-4), qt));
    if ~show_figure_flag, close(hfig); end
end