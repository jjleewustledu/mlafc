%% SL_fMRI_initialization
% 04/18/19
% KP

% Outputs from this code are: sl_fc_mean and sl_fc_gsp 

% sl_fc_mean: from 3. SL-fMRI mean map using 100 GSP subjects 
            % 2825 (# spheres) x 18611 (# grey matter voxels)            
% sl_fc_gsp : from 4. SL-fMRI correlation maps of 100 GSP subjects 
            % 100 (# subjects) x 65549 (# whole brain voxels)
            % this variable represents inter-subject variability across 
            % subjects in each voxel

            
% If you want to change the sphere radius or grid spacing, you should run
% this code again to create a new SL-fMRI "ground truth map"
% Section 2 (Store sphere voxels indices) - R (sphere radius); I (grid
% spacing); numvox (minimum number of voxels required to use the sphere)

%% 1. Initialization
clear; clc; close all

addpath(genpath('/data/nil-bluearc/shimony/Park/Matlab'))

% Load GLM mask, atlas image, and grey matter mask
[GMmask frames voxelsize] = read_4dfpimg_v1(['N21_aparc+aseg_GMctx_on_711-2V_333_avg_zlt0.5_gAAmask_v1.4dfp.img'],1,'littleendian');
[glmmask frames voxelsize] = read_4dfpimg_v1(['glm_atlas_mask_333.4dfp.img'],1,'littleendian');
[atlas frames voxelsize] = read_4dfpimg_v1('MNI152_T1_0.5mm_on_711-2B_333.4dfp.img',1,'littleendian');

% Configure masks and coordinate system
Nx = 48; Ny = 64; Nz = 48; % 333 space
N3D = Nx*Ny*Nz; 

glmmask(find(glmmask))=1; GMmsk_for_glm = GMmask(find(glmmask)); glmmask_3d = reshape(glmmask, [Nx, Ny, Nz]);
glmmsk333=reshape(glmmask,[48 64 48]);glmmsk333=flipdim(glmmsk333,2);GMmask333=reshape(GMmask,[48 64 48]);GMmask333=flipdim(GMmask333,2);
[ X1,X2,X3,xgrid ] = Generate_4dfpcoords( ); glmmaskind = find(glmmask);glmmaskind = find(glmmask);

%% 2. Store sphere voxels indices

R = 3; % sphere radius
I = R; % grid spacing
numvox = 10; % minimum # voxs required 

sphere_all = zeros(Nx, Ny, Nz);
ind = 1; 
for i = R:I:Nx-R
    for j = R:I:Ny-R
        for k = R:I:Nz-R

            xr = i; yr = j; zr = k;
            sphere_mask = zeros(Nx, Ny, Nz);
            
            for x = 1:Nx
                for y = 1:Ny
                    for z = 1: Nz
                        if ( (x - xr).^2 + (y - yr).^2 + (z - zr).^2 < R^2 )
                            sphere_mask(x,y,z) = 1; %//set elements within ellipsoid to 1
                        end
                    end
                end
            end
            
            sphere_mask = sphere_mask.*glmmask_3d;
            sphere_mask_glm = sphere_mask(glmmaskind);
            sphere_vox = find(sphere_mask_glm); % indices of sphere vox
            
            if length(sphere_vox) < numvox
                continue
            else
                sv{ind} = sphere_vox;
                ind = ind+1; 
            end
        end
    end
end

%% 3. SL-fMRI mean map using 100 GSP subjects 

% refdir = '/data/nil-bluearc/shimony/Park/GTM500Perceptron'; cd(refdir)
% refids = importdata('GTM500.lst');
% ext = '_faln_dbnd_xr3d_atl_g7_bpss_resid.mat';
% 
% sl_fcz_sum = zeros(2825,18611);
% for refnum = 1:100
% 
%     tic
%     refid = refids{refnum};
%     load(fullfile(refdir, refid, 'Perceptron', sprintf('%s%s', refid, ext)),'bold_frames')
%     
%     for ind = 1:length(sv)
%         
%         sphere_vox = sv{ind};
%         bf_sv = mean(bold_frames(:, sphere_vox),2);
%         sl_fc(ind, :) = corr_new(bf_sv, bold_frames(:, find(GMmsk_for_glm))); % SL_FC
%         
%     end
%     
%     % fisher-z-transform
%     sl_fcz = atanh(sl_fc);
%     sl_fcz_sum = sl_fcz_sum + sl_fcz; 
%         
%     fprintf('%s processed: %f\n', num2str(refnum), toc)
%     clear bold_frames sl_fc sl_fcz
%     
% end
% 
% % average
% sl_fcz_mean = sl_fcz_sum/100; 
% % tanh
% sl_fc_mean = tanh(sl_fcz_mean);

%% 4. SL-fMRI correlation maps of 100 GSP subjects 

% refdir = '/data/nil-bluearc/shimony/Park/GTM500Perceptron'; cd(refdir)
% refids = importdata('GTM500.lst');
% ext = '_faln_dbnd_xr3d_atl_g7_bpss_resid.mat';
% 
% for refnum = 1:100
% 
%     tic
%     refid = refids{refnum};
%     load(fullfile(refdir, refid, 'Perceptron', sprintf('%s%s', refid, ext)),'bold_frames')
%     
%     for ind = 1:length(sv)
%         
%         sphere_vox = sv{ind};
%         bf_sv = mean(bold_frames(:, sphere_vox),2);
%         sl_fc(ind, :) = corr_new(bf_sv, bold_frames(:, find(GMmsk_for_glm))); % SL_FC
%         
%     end
%     
%     r = corr_kp(sl_fc, sl_fc_mean);
% 
%     r_glm = zeros(1,65549);
%     count = zeros(1,65549);
% 
%     for ind = 1:length(sv)
% 
%         temp = zeros(1,65549); counttemp = zeros(1,65549);
%         temp(sv{ind}) = atanh(r(ind)); counttemp(sv{ind}) = 1;
%         r_glm = r_glm + temp; count = count + counttemp; 
% 
%     end
% 
%     sl_fc_gsp(refnum,:) = tanh(r_glm./count);
% 
%     fprintf('%s processed: %f\n', num2str(refnum), toc)
%     
% end