classdef AFC 
	%% AFC  
    % & \mathbb{E}_{\sigma' \in S, \tau' \in T} \qty\Bigg{
    %   \text{corr}_{\text{KP}} \qty\Bigg(
    %             \text{corr}\qty( \;
    %                        \mathbb{E}_{\qty|\vb*{\xi}^{\sigma'\tau'}|} \vb*{B}^{\tau'}(\vb*{\xi}^{\sigma'\tau'}), \;
    %                        \vb*{B}^{\tau'}(\vb*{\eta^{\tau'}}) ), \;
    %             \mathbb{E}_{\sigma \in S, \tau \in T} \qty\bigg{
    %             \text{corr}\qty\big( \;
    %                        \mathbb{E}_{\qty|\vb*{\xi}^{\sigma\tau}|}\vb*{B}^\tau(\vb*{\xi}^{\sigma\tau}), \;
    %                        \vb*{B}^\tau(\vb*{\eta}^\tau) ) }
    %             ) } \\ \\
    % & \text{for BOLD series $\vb*{B}^\tau$ of subject $\tau \in$ cohort $T$,} \\
    % & \text{voxels $\vb*{\xi}^{\sigma\tau}$ in sphere $\sigma \in$ searchlight $S$;} \\
    % & \text{voxels $\vb*{\eta}^\tau$ are samples of the subjects' cortex in atlas space}

	%  $Revision$
 	%  was created 29-May-2019 16:35:03 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlafc/src/+mlafc.
 	%% It was developed on Matlab 9.5.0.1067069 (R2018b) Update 4 for MACI64.  Copyright 2019 Kiyun Park, John Joowon Lee.
 	
	properties (Dependent)
        atlas_1d
 		registry
        sv
    end
    
    methods
        
        %% GET
        
        function g = get.atlas_1d(~)
            g = mlperceptron.PerceptronRegistry.read_MNI152_T1_0p5mm();
        end
        function g = get.registry(this)
            g = this.registry_;
        end
        function g = get.sv(this)
            if isempty(this.sv_)
                this = this.SL_fMRI_initialization();
            end
            g = this.sv_;
        end
        
        %%
        
        function MLPAFC(this)
            
            import mlpark.SearchLight
            import mlperceptron.Fourdfp
            
            addpath(genpath('Z:\shimony\Park\Matlab'))
            tgtatlas = 'Z:\shimony\Park\atlas\Targets\MNI152_711-2B_333.4dfp.img';
            tgtanatData = Fourdfp.Read4dfp(tgtatlas);
            
            %% MLP AFC
            numsub = 150;
            MLP_GTM_150 = this.registry.gtm_y; % GTM_y(:,:,1:numsub);
            
            ind = 1;
            for sub1 = 1:numsub
                for sub2 = sub1+1:numsub
                    
                    MLP_GTM_sub1 = MLP_GTM_150(:,:,sub1);
                    MLP_GTM_sub2 = MLP_GTM_150(:,:,sub2);
                    
                    % rmse
                    sq = (MLP_GTM_sub1 - MLP_GTM_sub2).^2;
                    m = sum(sq')/8;
                    r = sqrt(m);
                    mlp_rmse_con150(:,ind) = r; %#ok<AGROW>
                    
                    % cc
                    mlp_cc(:,ind) = SearchLight.corr_v2(MLP_GTM_sub1, MLP_GTM_sub2); %#ok<AGROW>
                    
                    ind = ind+1;
                end
            end
            
            %%
            E_mean = mean(mlp_rmse_con150');
            
            intersubject_var = zeros(1,147456);
            intersubject_var(this.glmmsk_indices_) = E_mean;
            
            figure; 
            SearchLight.MapOverlay(tgtanatData, intersubject_var,'caxvals',[0 0.2],'ttlvec', "Intersubject Variability: GSP 100")
            
            %%
            cc_mean = mean(mlp_cc');
            
            intersubject_var = zeros(1,147456);
            intersubject_var(this.glmmsk_indices_) = cc_mean;
            
            figure; 
            SearchLight.MapOverlay(tgtanatData, intersubject_var,'caxvals',[0 1],'ttlvec', "Intersubject Variability: GSP 100")
        end
        function tmap = mlp_afc(this, pat_mlp_v3)
            
            addpath(genpath('Z:\shimony\Park\Matlab'))
            load('MLPAFC/mlp_rmse_con100.mat')
            load('MLPAFC/MLP_GTM_100')
            
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
            
            for i = 1:length(this.GMmsk_for_glm_)
                
                [~,~,~,stats] = ttest2(mlp_rmse_con(i,:),mlp_rmse_pat(i,:));
                tmap(i) = stats.tstat;
                
            end
            
        end
        function SL_AFC(this, varargin)
            %% SL_AFC
            % 04/18/19
            % KP
            
            % Usage:
            % 0) Ensure completion of SL_fMRI_initialization
            % 1) Runs Multi-layer Perceptron (MLP) - Perceptron_Release
            % 2) Create a patient specific SL_fMRI map
            % 3) SL-derived aberrant functional connectivity            
            
            import mlperceptron.PerceptronRelease
            import mlpark.SearchLight
            import mlpark.*
            
            ip = inputParser;
            addRequired(ip, 'patientdir', @isfolder)
            addRequired(ip, 'patientid', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;   
            
            %% 0.  SL_fMRI_initialization            
            
            this.sv_ = this.sv;
            
            %% 1. MLP
            
            % Example:
            % patientdir = '\\bciserver2\Kay\Emily\New_Epilepsy\PT15';
            % patientid  = 'PT15';                     
            
            PerceptronRelease.Perceptron_Release(ipr.patientdir, ipr.patientid)
            
            %% 2. Create a patient specific SL-fMRI map
            
            %ext = '_faln_dbnd_xr3d_uwrp_atl_uout_resid.mat';
            %load (fullfile(patientdir, 'Perceptron', sprintf('%s%s', patientid, ext)), 'bold_frames')

            load(fullfile(ipr.patientdir, 'Perceptron', sprintf('%s%s', ipr.patientid, this.registry_.perceptron_resid_mat)), 'bold_frames')
            % bold_frames is 2307 x 65549 single, time-samples x parenchyma voxels; N_times is variable across studies
            
            sl_fc = zeros(length(this.sv), sum(this.GMmsk_for_glm_), class(bold_frames));
            for ind = 1:length(this.sv) % EXPENSIVE without prealloc; 1:2825
                fprintf('SL_AFC:  index %i of %i\n', ind, length(this.sv))
                sphere_vox = this.sv{ind}; % 13 x 1
                bf_sv = nanmean(bold_frames(:, sphere_vox), 2);
                sl_fc(ind, :) = PerceptronRelease.corr_new(bf_sv, bold_frames(:, find(this.GMmsk_for_glm_)));
                % sl_fc is 2825 x 18611, times-samples x found gray-matter voxels
            end
            
            r = SearchLight.corr_kp(sl_fc, this.sl_fc_mean_); % 1 x 2825 <= 2825 x 18611, 2825 x 18611
            
            r_glm = zeros(1, length(this.GMmsk_for_glm_)); % 1 x 65549
            count = zeros(1, length(this.GMmsk_for_glm_)); % 1 x 65549
            
            for ind = 1:length(this.sv)                
                temp = zeros(1, length(this.GMmsk_for_glm_)); 
                temp(this.sv{ind}) = atanh(r(ind)); 
                counttemp = zeros(1, length(this.GMmsk_for_glm_));
                counttemp(this.sv{ind}) = 1;
                r_glm = r_glm + temp; 
                count = count + counttemp;                
            end
            
            sl_fmri_pat = tanh(r_glm./count); % 1 x 65549
            
            %% 3. SL-AFC            
            
            afc_map = threshed_afc_map(sl_fmri_pat);
            
            %% Visualization
            
            figure
            SearchLight.MapOverlay(this.atlas_1d', afc_map)
            
            %% save afc_map file
            
            workdir = fullfile(ipr.patientdir, ['AFCmap_' datestr(now,30)]);
            mkdir(workdir)
            pwd0 = pushd(workdir);
            save('afc_map', 'afc_map')
            saveFigures
            popd(pwd0)
            
            function map = threshed_afc_map(sl_fmri_pat)                
                sd_thr = 3; % z-thresh for visualization   
                
                sl_fc_gsp_mean = nanmean(this.sl_fc_gsp_);              % 1 x 65549
                sl_fc_gsp_std  = nanstd(this.sl_fc_gsp_);               % 1 x 65549
                sl_fc_gsp_thr  = sl_fc_gsp_mean - sd_thr*sl_fc_gsp_std; % 1 x 65549

                afc_vox             = find(sl_fmri_pat < sl_fc_gsp_thr);     % 1 x 56281
                afc_glmmap          = zeros(1, length(this.GMmsk_for_glm_)); % 1 x 65549
                afc_glmmap(afc_vox) = 1;                                     % 1 x 65549            
                [~,~,~,N3D] = this.registry_.atlas_dims();
                map = zeros(1, N3D); 
                map(this.glmmsk_indices_) = afc_glmmap; % 1 x 147456
            end
        end
        function this = SL_fMRI_initialization(this, varargin)
            %% SL_fMRI_initialization
            %  @param init is logical
            %  @param make_sv is logical
            %  @param make_sl_fc_mean is logical
            %  @param make_sl_fc_gsp is logical
            %  @param sphere_radius
            %  @param grid_spacing
            %  @param min_num_vox
            
            % 04/18/19 KP
            % 06/06/19 JJL
            %
            % Outputs from this code are: this.sv_, this.sl_fc_mean_ and this.sl_fc_gsp_
            %
            % this.sl_fc_mean_: from 3. SL-fMRI mean map using 100 GSP subjects
            %             2825 (# spheres) x 18611 (# grey matter voxels)
            % this.sl_fc_gsp_ : from 4. SL-fMRI correlation maps of 100 GSP subjects
            %             100 (# subjects) x 65549 (# whole brain voxels)
            %             this variable represents inter-subject variability across subjects in each voxel            
            %
            % If you want to change the sphere radius or grid spacing, you should run
            % this code again to create a new SL-fMRI "ground truth map"
            
            % addpath(genpath('/data/nil-bluearc/shimony/Park/Matlab'))
            
            import mlperceptron.*
            
            ip = inputParser;
            addParameter(ip, 'init', true, @islogical)
            addParameter(ip, 'make_sv', true, @islogical)
            addParameter(ip, 'make_sl_fc_mean', false, @islogical)
            addParameter(ip, 'make_sl_fc_gsp', false, @islogical)
            addParameter(ip, 'sphere_radius', 3, @isnumeric)
            addParameter(ip, 'grid_spacing', 3, @isnumeric)
            addParameter(ip, 'min_num_vox', 10, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.sv_ = {};
            this = mlafc.AFC();
            
            %% 1. Initialization
            
            %% 2. Store sphere voxels indices
            %  numvox (minimum number of voxels required to use the sphere)

            if ipr.make_sv
                [Nx,Ny,Nz] = this.registry_.atlas_dims;
                glmmsk_3d = reshape(this.glmatl_, [Nx, Ny, Nz]);
                R = ipr.sphere_radius;
                I = ipr.grid_spacing;
                ind = 1;
                for xr = R:I:Nx-R
                    for yr = R:I:Ny-R
                        for zr = R:I:Nz-R

                            sphere_mask = zeros(Nx, Ny, Nz);
                            for x = 1:Nx
                                for y = 1:Ny
                                    for z = 1: Nz
                                        if ( (x - xr)^2 + (y - yr)^2 + (z - zr)^2 < R^2 )
                                            sphere_mask(x,y,z) = 1; % set elements within ellipsoid to 1
                                        end
                                    end
                                end
                            end

                            sphere_mask = sphere_mask .* glmmsk_3d;
                            sphere_mask_glm = sphere_mask(this.glmmsk_indices_);
                            sphere_vox = find(sphere_mask_glm); % indices => sphere_vox
                            if length(sphere_vox) >= ipr.min_num_vox
                                this.sv_{ind} = sphere_vox; 
                                ind = ind+1;
                            end
                        end
                    end
                end
            end
            
            %% 3. SL-fMRI mean map using 100 GSP subjects
            
            if ipr.make_sl_fc_mean
                import mlperceptron.PerceptronRelease;
                pwd0 = pushd(this.registry_.gtm500_dir);           
                found_GMmsk_for_glm = find(this.GMmsk_for_glm_); % 18611 x 1
                sl_fcz_sum = zeros(length(this.sv_), length(found_GMmsk_for_glm)); 
                
                for refnum = 1:this.registry_.ref_count
                    
                    tic
                    refid = this.registry_.gtm500_ids{refnum};
                    load(fullfile(this.registry_.gtm500_dir, refid, 'Perceptron', sprintf('%s%s', refid, this.registry_.ref_resid_mat)),'bold_frames')
                    % bold_frames has 245 x 65549
                    for ind = 1:length(this.sv_)                        
                        sphere_vox = this.sv_{ind}; % 13 x 1
                        bf_sv = mean(bold_frames(:, sphere_vox),2); % 245 x 1
                        sl_fc(ind, :) = PerceptronRelease.corr_new(bf_sv, bold_frames(:, found_GMmsk_for_glm)); %#ok<AGROW> % SL_FC
                        % 2825 x 18611
                    end
                    
                    % Fisher z-transform
                    sl_fcz = atanh(sl_fc);
                    sl_fcz_sum = sl_fcz_sum + sl_fcz;
                    
                    fprintf('%s processed: %f\n', num2str(refnum), toc)
                    clear bold_frames sl_fc sl_fcz                    
                end
                
                % average
                sl_fcz_mean = sl_fcz_sum / this.registry_.ref_count;
                % tanh <-> inverse Fisher z-transform
                sl_fc_mean = tanh(sl_fcz_mean);
                save(this.registry_.sl_fc_mean_mat, 'sl_fc_mean')
                this.sl_fc_mean_ = sl_fc_mean;
                clear('sl_fc_mean')
                popd(pwd0)
            end
            
            %% 4. SL-fMRI correlation maps of 100 GSP subjects
            
            if ipr.make_sl_fc_gsp
                import mlperceptron.PerceptronRelease;
                import mlpark.SearchLight;
                pwd0 = pushd(this.registry_.gtm500_dir);                 
                sl_fc_gsp = zeros(this.registry_.ref_count, length(this.GMmsk_for_glm_));
                
                for refnum = 1:this.registry_.ref_count
                    
                    tic
                    refid = this.registry_.gtm500_ids{refnum};
                    load(fullfile(this.registry_.gtm500_dir, refid, 'Perceptron', sprintf('%s%s', refid, this.registry_.ref_resid_mat)),'bold_frames')
                    
                    for ind = 1:length(this.sv_)
                        
                        sphere_vox = this.sv_{ind};
                        bf_sv = mean(bold_frames(:, sphere_vox),2);
                        sl_fc(ind, :) = PerceptronRelease.corr_new(bf_sv, bold_frames(:, find(this.GMmsk_for_glm_))); % SL_FC
                        
                    end
                    
                    r = SearchLight.corr_kp(sl_fc, this.sl_fc_mean_);
                    
                    r_glm = zeros(1, length(this.GMmsk_for_glm_));
                    count = zeros(1, length(this.GMmsk_for_glm_));                    
                    for ind = 1:length(this.sv_)                        
                        temp = zeros(1, length(this.GMmsk_for_glm_)); 
                        counttemp = zeros(1, length(this.GMmsk_for_glm_));
                        temp(this.sv_{ind}) = atanh(r(ind)); 
                        counttemp(this.sv_{ind}) = 1;
                        r_glm = r_glm + temp; 
                        count = count + counttemp;                        
                    end
                    
                    % tanh <-> inverse Fisher z-transform
                    sl_fc_gsp(refnum,:) = tanh(r_glm./count);
                    fprintf('%s processed: %f\n', num2str(refnum), toc)                    
                end
                save(this.registry_.sl_fc_gsp_mat, 'sl_fc_gsp')
                this.sl_fc_gsp_ = sl_fc_gsp;
                clear('sl_fc_gsp')
                popd(pwd0)
            end
        end
        
        function sl_fc_gsp = expectation_corr_of_corr_bold(this)
            
            
        end
		  
 		function this = AFC(varargin)
 			%% AFC
 			%  @param .

            import mlperceptron.*;
            
            this.registry_ = mlafc.AFCRegistry.instance();            
            if isfile(this.registry_.sl_fc_mean_mat)
                load(this.registry_.sl_fc_mean_mat)
                this.sl_fc_mean_ = sl_fc_mean;
                clear('sl_fc_mean')
            end
            if isfile(this.registry_.sl_fc_gsp_mat)
                load(this.registry_.sl_fc_gsp_mat)
                this.sl_fc_gsp_ = sl_fc_gsp;
                clear('sl_fc_gsp')
            end
            
            this.glmatl_ = PerceptronRegistry.read_glm_atlas_mask();
            this.glmatl_(find(this.glmatl_)) = 1;  
            this.glmmsk_indices_ = find(this.glmatl_);
            GMctx = PerceptronRegistry.read_N21_aparc_aseg_GMctx(); % 147456 x 1
            this.GMmsk_for_glm_ = GMctx(find(this.glmatl_)); %#ok<*FNDSB> 
 		end
    end 
    
    %% HIDDEN
    
    properties (Hidden)
        glmatl_         % 147456 x 1
        glmmsk_indices_ % 65549  x 1
        GMmsk_for_glm_  % 65549  x 1
        registry_
        sl_fc_mean_     % 2825 x 18611
        sl_fc_gsp_      % 100  x 65549
        sv_             % 1    x 2825  cells with variable size [(x > 1) 1]
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

