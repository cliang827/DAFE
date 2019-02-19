function auc_result = cmc2auc(cmc_result)


auc_result = 0.5*(2*sum(cmc_result) - cmc_result(1) - cmc_result(end))/(length(cmc_result)-1);