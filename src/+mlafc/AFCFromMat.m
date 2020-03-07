classdef AFCFromMat < mlafc.AFC
	%% AFCFROMMAT  

	%  $Revision$
 	%  was created 17-Jan-2020 01:00:53 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlafc/src/+mlafc.
 	%% It was developed on Matlab 9.7.0.1261785 (R2019b) Update 3 for MACI64.  Copyright 2020 John Joowon Lee.
 	
	properties (Constant)
 		MSC_SUBJECTS_MAT = '/data/nil-bluearc/shimony/jjlee/MSC_Avis_denoising/MSC_subjects.mat'
 	end

	methods 
        function this = makeMscMap(this, varargin)
            %% MAKEMSCMAP
            %  Preconditions:
            %  -- completion of SL_fMRI_initialization
            %  -- completion of mlperceptron.PerceptronRelease factory
            %  @param patientdir is a folder.
            %  @param patientid is char.
            %  @param afc_filename is char.
            %  @param Nframes is numeric.
            %  @param feature is a filename, NIfTI or 4dfp.
            %  @returns instance of mlafc.AFC.  
            %  @returns inputParser.Results.            
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'patientdir', @isfolder)
            addRequired( ip, 'patientid', @ischar)
            addParameter(ip, 'afc_filename', ['AFC_this_' datestr(now, 30) '.mat'], @ischar)
            addParameter(ip, 'feature_filename', '', @isfile)
            addParameter(ip, 'load_afc', true, @islogical)
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.patientdir = ipr.patientdir;
            this.patientid = ipr.patientid;
            
            if ipr.load_afc && isfile(ipr.afc_filename)
                load(ipr.afc_filename, 'this')
            else
                this = this.makeSearchlightMap(varargin{:});
            end
            
            workdir = fileparts(ipr.afc_filename);
            mkdir(workdir)
            pwd0 = pushd(workdir);
            save(this.product)
            save(ipr.afc_filename, 'this',  'ipr') 
            figure
            this.feature = zeros(1, 147456);
            this.mapOverlay2(this.read_T1_333', this.afc_map, this.feature)                
            saveFigures(workdir)
            popd(pwd0)
        end
        function [this,ipr] = makeSearchlightMap(this, varargin)
            %% MAKESEARCHLIGHTMAP
            %  Preconditions:
            %  -- completion of SL_fMRI_initialization
            %  -- completion of mlperceptron.PerceptronRelease factory
            %  @param patientdir is a folder.
            %  @param patientid is char.
            %  @param afc_filename is char.
            %  @param Nframes is numeric.
            %  @returns instance of mlafc.AFC.
            %  @returns inputParser.Results.
            
            import mlperceptron.PerceptronFromMat
            import mlpark.SearchLight
            import mlpark.*
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'patientdir', @isfolder)
            addRequired( ip, 'patientid', @ischar)
            addParameter(ip, 'afc_filename', ['AFC_this_' datestr(now, 30) '.mat'], @ischar)
            addParameter(ip, 'Nframes', [], @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;             
            disp(ipr)
            this.patientdir = ipr.patientdir;
            this.patientid = ipr.patientid;
            
            %% using results of SL_fMRI_initialization            
            
            this.sv_ = this.sv;
            
            %% using results of mlperceptron.PerceptronRelease factory                   
            
            PerceptronFromMat.createMscProbMaps(6, 'mscSubjectsMat', this.MSC_SUBJECTS_MAT)
            
            %% 

            load(fullfile(this.patientdir, 'Perceptron', sprintf('%s%s', this.patientid, this.registry.perceptron_resid_mat)), 'bold_frames')
            % bold_frames is 2307 x 65549 single, time-samples x parenchyma voxels; N_times is variable across studies
            
            if ~isempty(ipr.Nframes)
                bold_frames = bold_frames(1:ipr.Nframes, :);
            end 
            fprintf('makeSearchlightMap.bold_frames.size -> %s\n', mat2str(size(bold_frames)))
            
            sl_fc = zeros(length(this.sv), sum(this.GMmsk_for_glm_), class(bold_frames));
            for ind = 1:length(this.sv) % EXPENSIVE without prealloc; 1:2825
                sphere_vox = this.sv{ind}; % 13 x 1
                bf_sv = nanmean(bold_frames(:, sphere_vox), 2);
                sl_fc(ind, :) = this.similarity(bf_sv, bold_frames(:, find(this.GMmsk_for_glm_)), 'kind', 'ch');
                % sl_fc is 2825 x 18611, times-samples x found gray-matter voxels
                
                if 0 == mod(ind, 100)
                    fprintf('makeSearchlightMap:  assigned sphere_vox for index %i of %i\n', ind, length(this.sv))
                end
            end
            
            r = this.similarity(sl_fc, this.sl_fc_mean_); % 1 x 2825 <= 2825 x 18611, 2825 x 18611
            
            r_glm = zeros(1, length(this.GMmsk_for_glm_)); % 1 x 65549
            count = zeros(1, length(this.GMmsk_for_glm_)); % 1 x 65549            
            for ind = 1:length(this.sv)                
                temp = zeros(1, length(this.GMmsk_for_glm_)); 
                temp(this.sv{ind}) = atanh(r(ind)); 
                counttemp = zeros(1, length(this.GMmsk_for_glm_));
                counttemp(this.sv{ind}) = 1;
                r_glm = r_glm + temp; 
                count = count + counttemp;
                
                if 0 == mod(ind, 100)
                    fprintf('makeSearchlightMap:  assigned count for index %i of %i\n', ind, length(this.sv))
                end                
            end            
            this.sl_fmri_pat = tanh(r_glm./count); % 1 x 65549
            this.afc_map = this.diff_map(this.sl_fmri_pat);
        end
        function [img,frames,voxelsize] = read_T1_333(this)
            import mlperceptron.*
            ss = strsplit(this.patientid, '_subject');
            numid = [ss{1} sprintf('%02i', str2double(ss{2}))];
            [img,frames,voxelsize] = Fourdfp.read_4dfpimg_v1( ...
                fullfile(this.patientdir, 'atlas', [numid '_mpr_avgT_333_t88.4dfp.img']));
        end
		  
 		function this = AFCFromMat(varargin)
 			%% AFCFROMMAT
 			%  @param .

 			this = this@mlafc.AFC(varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

