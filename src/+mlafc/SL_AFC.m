%% SL_AFC
% 04/18/19
% KP

% Usage: 
% 1) Runs Multi-layer Perceptron (MLP) - Perceptron_Release
% 2) Create a patient specific SL_fMRI map
% 3) SL-derived aberrant functional connectivity

%% 1. MLP

% Example:
% patientdir = '\\bciserver2\Kay\Emily\New_Epilepsy\PT15';
% patientid  = 'PT15';

Perceptron_Release(patientdir,patientid)

%% 2. Create a patient specific SL-fMRI map

load('sl_fc_mean.mat')

ext = '_faln_dbnd_xr3d_uwrp_atl_uout_resid.mat';
load (fullfile(patientdir, 'Perceptron', sprintf('%s%s', patientid, ext)), 'bold_frames')
        
for ind = 1:length(sv)

    sphere_vox = sv{ind};
    bf_sv = nanmean(bold_frames(:, sphere_vox),2);
    sl_fc(ind, :) = corr_new(bf_sv, bold_frames(:, find(GMmsk_for_glm))); % SL_FC

end

r = corr_kp(sl_fc, sl_fc_mean);

r_glm = zeros(1,65549);
count = zeros(1,65549);

for ind = 1:length(sv)

    temp = zeros(1,65549); counttemp = zeros(1,65549);
    temp(sv{ind}) = atanh(r(ind)); counttemp(sv{ind}) = 1;
    r_glm = r_glm + temp; count = count + counttemp; 

end

sl_fmri_pat = tanh(r_glm./count); % size: 1 x 65549

%% 3. SL-AFC

load('sl_fc_gsp.mat')

sd_thr = 3; 
sl_fc_gsp_mean = nanmean(sl_fc_gsp);
sl_fc_gsp_std = nanstd(sl_fc_gsp);
sl_fc_gsp_thr = sl_fc_gsp_mean - sd_thr*sl_fc_gsp_std;

afc_glmmap = zeros(1,65549);
afc_vox = find(sl_fmri_pat < sl_fc_gsp_thr);
afc_glmmap(afc_vox) = 1;

afc_map = zeros(1,N3D); afc_map(glmmaskind) = afc_glmmap; 

%% Visualization

figure; MapOverlay(atlas, afc_map)


%% save afc_map file

mkdir(sprintf('%s\\AFCmap', patientdir))
cd(sprintf('%s\\AFCmap', patientdir))
save('afc_map', 'afc_map')

