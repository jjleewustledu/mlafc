function tmap = mlp_afc(pat_mlp_v3)

addpath(genpath('Z:\shimony\Park\Matlab'))
load('mlp_rmse_con100')
load('MLP_GTM_100')

MLP_GTM = MLP_GTM_100;
mlp_rmse_con = mlp_rmse_con100;

pat_mlp_y = 1./(1 + exp(-pat_mlp_v3));
pat_mlp = pat_mlp_y./sum(pat_mlp_y')';

for sub = 1:50
    MLP_GTM_sub = MLP_GTM(:,:,sub);
    sq = (pat_mlp - MLP_GTM_sub).^2;
    m = (sum(sq')/8);
    r = sqrt(m);
    mlp_rmse_pat(:,sub) = r; 
end

for i = 1:65549
    
    [~,~,~,stats] = ttest2(mlp_rmse_con(i,:),mlp_rmse_pat(i,:));
    tmap(i) = stats.tstat;
    
end

end
