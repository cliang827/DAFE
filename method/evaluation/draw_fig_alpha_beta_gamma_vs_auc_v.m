%% figure 2: alpha/beta/gamma vs. auc(f) and mean(v)/std(v)
hfig = figure;
assert(ctrl_para.exp_para.tot_query_times==3);

alpha_auc = alpha_set;
alpha_v_mean = zeros(size(alpha_set,1),3);
alpha_v_std = zeros(size(alpha_set,1),3);
for i=1:length(alpha_set)
    alpha_auc(i) = mean(auc_f_kmean(alpha==alpha_set(i),2));
    
    alpha_v_mean(i,:) = mean(v_max(alpha==alpha_set(i),:));
    alpha_v_std(i,:) = mean(v_std(alpha==alpha_set(i),:));
end
subplot(3,3,1);
plot(alpha_set, alpha_auc, 'r-*'); grid on;
xlabel('log(\alpha)'); ylabel('auc'); title('\alpha-auc');

subplot(3,3,2);
plot(alpha_set, alpha_v_mean(:,1), 'g--o'); hold on;
plot(alpha_set, alpha_v_mean(:,2), 'b-.s'); hold on; 
plot(alpha_set, alpha_v_mean(:,3), 'r-*'); hold on; grid on;
xlabel('log(\alpha)'); ylabel('mean(v)'); title('\alpha-v-mean');
legend({'qt=1','qt=2','qt=3'});

subplot(3,3,3);
plot(alpha_set, alpha_v_std(:,1), 'g--o'); hold on;
plot(alpha_set, alpha_v_std(:,2), 'b-.s'); hold on; 
plot(alpha_set, alpha_v_std(:,3), 'r-*'); hold on; grid on;
xlabel('log(\alpha)'); ylabel('std(v)'); title('\alpha-v-std');
legend({'qt=1','qt=2','qt=3'});


beta_auc = beta_set;
beta_v_mean = zeros(size(beta_set,1),3);
beta_v_std = zeros(size(beta_set,1),3);
for i=1:length(beta_set)
    beta_auc(i) = mean(auc_f_kmean(beta==beta_set(i),2));
    
    beta_v_mean(i,:) = mean(v_max(beta==beta_set(i),:));
    beta_v_std(i,:) = mean(v_std(beta==beta_set(i),:));
end
subplot(3,3,4);
plot(beta_set, beta_auc, 'r-*'); 
xlabel('\beta%'); ylabel('auc'); title('\beta%-auc'); grid on;

subplot(3,3,5);
plot(beta_set, beta_v_mean(:,1), 'g--o'); hold on;
plot(beta_set, beta_v_mean(:,2), 'b-.s'); hold on; 
plot(beta_set, beta_v_mean(:,3), 'r-*'); hold on;  grid on;
xlabel('\beta%'); ylabel('mean(v)'); title('\beta%-v-mean'); grid on;
legend({'qt=1','qt=2','qt=3'});

subplot(3,3,6);
plot(beta_set, beta_v_std(:,1), 'g--o'); hold on;
plot(beta_set, beta_v_std(:,2), 'b-.s'); hold on; 
plot(beta_set, beta_v_std(:,3), 'r-*'); hold on;  grid on;
xlabel('\beta%'); ylabel('std(v)'); title('\beta%-v-std'); grid on;
legend({'qt=1','qt=2','qt=3'});


gamma_auc = gamma_set;
gamma_v_mean = zeros(size(gamma_set,1),3);
gamma_v_std = zeros(size(gamma_set,1),3);
for i=1:length(gamma_set)
    gamma_auc(i) = mean(auc_f_kmean(gamma==gamma_set(i),2));
    
    gamma_v_mean(i,:) = mean(v_max(gamma==gamma_set(i),:));
    gamma_v_std(i,:) = mean(v_std(gamma==gamma_set(i),:));
end
subplot(3,3,7);
plot(gamma_set, gamma_auc, 'r-*'); 
xlabel('log(\gamma)'); ylabel('auc'); title('\gamma-auc'); grid on;

subplot(3,3,8);
plot(gamma_set, gamma_v_mean(:,1), 'g--o'); hold on; 
plot(gamma_set, gamma_v_mean(:,2), 'b-.s'); hold on; 
plot(gamma_set, gamma_v_mean(:,3), 'r-*'); hold on; grid on;
xlabel('log(\gamma)'); ylabel('mean(v)'); title('\gamma-v-mean'); grid on;
legend({'qt=1','qt=2','qt=3'});

subplot(3,3,9);
plot(gamma_set, gamma_v_std(:,1), 'g--o'); hold on;
plot(gamma_set, gamma_v_std(:,2), 'b-.s'); hold on;
plot(gamma_set, gamma_v_std(:,3), 'r-*'); hold on; grid on;
xlabel('log(\gamma)'); ylabel('std(v)'); title('\gamma-v-std'); grid on;
legend({'qt=1','qt=2','qt=3'});

saveas(hfig, [result_mat_file(1:end-4), '-fV.fig']);
if ~show_figure_flag, close(hfig); end