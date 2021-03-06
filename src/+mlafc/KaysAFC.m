classdef KaysAFC < AFC
	%% KAYSAFC  

	%  $Revision$
 	%  was created 01-May-2021 18:22:35 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlafc/src/+mlafc.
 	%% It was developed on Matlab 9.10.0.1602886 (R2021a) for MACI64.  Copyright 2021 John Joowon Lee.
 	
	properties (Dependent)
        atlas_1d
    end
    
	properties 	
        feature	
        similarityKind = 'kp'
        sl_fc_gsp_ % 100  x 65549 ~ |GSP| x |parenchyma|
        sl_fmri_pat
    end
    
    methods (Static)
        function [sim,featifc,funcifc] = calcdice(featdata, funcdata, varargin)
            import mlafc.AFC.removeInfratentorial
            
            ip = inputParser;
            addRequired(ip, 'featdata', @isnumeric)
            addRequired(ip, 'funcdata', @isnumeric)   
            addParameter(ip, 'Nsigma', 3, @isnumeric)
            parse(ip, featdata, funcdata, varargin{:})
            ipr = ip.Results;
            
            [GLMmask,~,~,glmic] = mlperceptron.PerceptronRegistry.read_glm_atlas_mask();
            GLMmask(find(GLMmask)) = 1; 
            GLMmask = removeInfratentorial(GLMmask);
            
            % reshape
            GLMmask  = reshape(GLMmask,  [48 64 48]);
            featdata = reshape(featdata, [48 64 48]) .* GLMmask; 
            funcdata = reshape(funcdata, [48 64 48]) .* GLMmask;
            
            % threshold funcdata
            mu = nanmean(funcdata(logical(GLMmask)));
            sigma = nanstd(funcdata(logical(GLMmask)));
            funcdata = funcdata < mu - ipr.Nsigma*sigma;
            
            featifc = copy(glmic.fourdfp);
            featifc.img = single(featdata);
            featifc.filepath = pwd;
            featifc.fileprefix = 'calcdice_label';
            funcifc = copy(glmic.fourdfp);
            funcifc.img = single(funcdata);
            funcifc.filepath = pwd;
            funcifc.fileprefix = 'calcdice_afc_3sd';
            
            sim = dice(logical(featdata), logical(funcdata));  
            fprintf('mlafc.AFC.calcdice.sim -> %g\n', sim)
        end
        function fsleyes(img, varargin)
            ip = inputParser;
            addParameter(ip, 'changeSign', false) % for viz.
            addParameter(ip, 'nanToZero', false)
            addParameter(ip, 'filename', '') % save on the fly while debugging
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            if ipr.changeSign
                img = -img;
            end
            if ipr.nanToZero
                img(isnan(img)) = 0;
            end
            if ismatrix(img)
                img = reshape(img, [48 64 48]);
            end
            img = flip(flip(img, 1), 2);
            ic  = mlfourd.ImagingContext2(img);
            ic.fsleyes
            if ~isempty(ipr.filename)
                ic.filename = ipr.filename;
                ic.save
            end
        end
        function d = kldiv(x, y, varargin)
            %% Kullback-Leibler or Jensen-Shannon divergence implemented by David Fass, 2016.
            %  @param required x, y are numeric.
            %  @param optional 'js' requests the J-S divergence.
            %  @param optional 'sym' requests a symmetrized K-L
            
            p = preprocess(x);
            q = preprocess(y);
            d = kldiv((1:length(p))', p, q, varargin{:});
            
            function p = preprocess(x)
                %% p := x interpreted as vector of probabilities
                
                if ndims(x) > 1
                    x = reshape(x, [numel(x) 1]);
                end
                x = ascol(x);
                x = x - min(x);
                p = x / sum(x) + eps;
            end
        end
        function mapOverlay(anatdata, funcdata, varargin)
            %% Adapted from codes by Kay Park
            
            import mlafc.AFC
            import mlpark.SearchLight
            import mlperceptron.Fourdfp
            import mlperceptron.PerceptronRegistry
            
            ip = inputParser;
            addRequired( ip, 'anatdata', @(x) true);
            addRequired( ip, 'funcdata', @(x) true);            
            addParameter(ip, 'ttlvec', ' ' , @(s)isstring(s));
            addParameter(ip, 'caxvals', [-.25, .25], @isnumeric)
            addParameter(ip, 'cmapopt', 'parula', @isstring)
            addParameter(ip, 'spacedim', '333', @isstring)
            addParameter(ip, 'plane', 'axial', @isstring)
            addParameter(ip, 'layoutsize', [4, 12], @isnumeric)
            addParameter(ip, 'snapsize', [4, 12], @isnumeric)
            addParameter(ip, 'snapind', 1:48, @isnumeric)
            parse(ip, anatdata, funcdata, varargin{:});
            ipr = ip.Results;
            
            % GLMmask
            GLMmask = PerceptronRegistry.read_glm_atlas_mask();
            GLMmask(find(GLMmask)) = 1; 
            
            % reshape data
            N3D = mlafc.AFCRegistry.instance().atlas_numel();
            anatdata = reshape(anatdata, [1, N3D]); % 1 x 147456 single
            funcdata = reshape(funcdata, [1, N3D]); 
            GLMmask  = reshape(GLMmask,  [1, N3D]); 
            
            % anatdata, funcdata
            anat_image_temp = SearchLight.rotflrst_image(anatdata, ipr.layoutsize, ipr.plane);     % 256 x 576  single
            funcdata_glm = funcdata .* GLMmask;                                              % 1 x 147456 single
            func_image_temp = SearchLight.rotflrst_image(funcdata_glm, ipr.layoutsize, ipr.plane); % 256 x 576  single
            
            if ~isequal(ipr.snapsize, ipr.layoutsize)
                anat_image = AFC.snapImage(anat_image_temp, ipr.snapsize, ipr.snapind, ipr.plane);
                func_image = AFC.snapImage(func_image_temp, ipr.snapsize, ipr.snapind, ipr.plane);
            else
                anat_image = anat_image_temp; % BUG:  noise [-realmax realmax]
                func_image = func_image_temp; % BUG:  noise [0 1]
            end
            
            pboptions = [size(anat_image, 2), size(anat_image, 1),1];
            alphaset = 0.8*ones(1, length(anat_image(:))); 
            alphaset(find(func_image(:) == 0)) = 0;           % 1 x 147456
            alphaset3d = reshape(alphaset, size(func_image)); % 256 x 576, [0 1]
            
            % plot
            ax1 = axes; 
            imagesc(anat_image); 
            pbaspect(pboptions); 
            axis off;
            % caxis([500 max(anat_image(:))-500])
            % caxis([0 max(anat_image(:))-150])
            title(ipr.ttlvec, 'interpreter','none')
            ax2 = axes; 
            imagesc(func_image); 
            alpha(alphaset3d); 
            pbaspect(pboptions);
            linkaxes([ax1, ax2])
            ax2.Visible = 'Off';
            colormap(ax1, 'gray');
            colormap(ax2, ipr.cmapopt)
            caxis(ipr.caxvals)
        end
        function mapOverlay2(anatdata, funcdata, featdata, varargin)
            %% Adapted from codes by Kay Park
            
            import mlafc.AFC
            import mlpark.SearchLight
            import mlperceptron.Fourdfp
            import mlperceptron.PerceptronRegistry
            
            ip = inputParser;
            addRequired( ip, 'anatdata', @(x) true)
            addRequired( ip, 'funcdata', @(x) true)
            addRequired( ip, 'featdata', @(x) true)
            addParameter(ip, 'ttlvec', ' ' , @(s)isstring(s))
            addParameter(ip, 'caxvals', [-.25, .25], @isnumeric)
            addParameter(ip, 'cmapopt', 'parula', @isstring)
            addParameter(ip, 'spacedim', '333', @isstring)
            addParameter(ip, 'plane', 'axial', @isstring)
            addParameter(ip, 'layoutsize', [4, 12], @isnumeric)
            addParameter(ip, 'snapsize', [4, 12], @isnumeric)
            addParameter(ip, 'snapind', 1:48, @isnumeric)
            parse(ip, anatdata, funcdata, featdata, varargin{:})
            ipr = ip.Results;
            
            % GLMmask
            GLMmask = PerceptronRegistry.read_glm_atlas_mask();
            GLMmask(find(GLMmask)) = 1; 
            
            % reshape data
            N3D = mlafc.AFCRegistry.instance().atlas_numel();
            GLMmask  = reshape(GLMmask,  [1, N3D]); % 1 x 147456 single 
            anatdata = reshape(anatdata, [1, N3D]);
            funcdata = reshape(funcdata, [1, N3D]); 
            featdata = reshape(featdata, [1, N3D]); 
            
            % anatdata, funcdata, featdata
            funcdata_glm = funcdata .* GLMmask;                                                    % 1 x 147456 single            
            featdata_glm = featdata .* GLMmask;   
            anat_image_temp = SearchLight.rotflrst_image(anatdata, ipr.layoutsize, ipr.plane);     % 256 x 576  single
            func_image_temp = SearchLight.rotflrst_image(funcdata_glm, ipr.layoutsize, ipr.plane);
            feat_image_temp = edge( ...
                              SearchLight.rotflrst_image(featdata_glm, ipr.layoutsize, ipr.plane));
            
            if ~isequal(ipr.snapsize, ipr.layoutsize)
                anat_image = AFC.snapImage(anat_image_temp, ipr.snapsize, ipr.snapind, ipr.plane);
                func_image = AFC.snapImage(func_image_temp, ipr.snapsize, ipr.snapind, ipr.plane);
                feat_image = AFC.snapImage(feat_image_temp, ipr.snapsize, ipr.snapind, ipr.plane);
            else
                anat_image = anat_image_temp; % BUG:  noise [-realmax realmax]
                func_image = func_image_temp; % BUG:  noise [0 1]
                feat_image = feat_image_temp; % BUG:  noise [0 1]
            end
            
            pboptions = [size(anat_image, 2), size(anat_image, 1), 1];
            alphaFeat = ones(1, length(anat_image(:))); 
            alphaFeat(find(feat_image(:) == 0)) = 0;           % 1 x 147456
            alphaFeat = reshape(alphaFeat, size(func_image)); % 256 x 576, [0 1]
            alphaFunc = ones(1, length(anat_image(:))); 
            alphaFunc(find(func_image(:) == 0)) = 0;           % 1 x 147456
            alphaFunc = reshape(alphaFunc, size(func_image)); % 256 x 576, [0 1]
            
            % plot
            ax1 = axes; 
            imagesc(anat_image); 
            pbaspect(pboptions); 
            colormap(ax1, 'gray');
            axis off;
            title(ipr.ttlvec, 'interpreter','none')
            
            ax3 = axes; 
            imagesc(func_image); 
            alpha(alphaFunc); 
            pbaspect(pboptions);
            linkaxes([ax1, ax3])
            ax3.Visible = 'Off';
            colormap(ax3, ipr.cmapopt)
            caxis(ipr.caxvals)
            
            ax2 = axes; 
            imagesc(feat_image); 
            alpha(alphaFeat); 
            pbaspect(pboptions);
            linkaxes([ax1, ax2])
            ax2.Visible = 'Off';
            colormap(ax2, 'gray');
            caxis(ipr.caxvals)
        end
        function plotroc(featdata, funcdata, varargin)
            import mlafc.AFC.removeInfratentorial
            
            ip = inputParser;
            addRequired(ip, 'featdata', @isnumeric)
            addRequired(ip, 'funcdata', @isnumeric)            
            parse(ip, featdata, funcdata, varargin{:})
            GLMmask = mlperceptron.PerceptronRegistry.read_glm_atlas_mask();
            GLMmask(find(GLMmask)) = 1; 
            GLMmask = removeInfratentorial(GLMmask);
            N3D = mlafc.AFCRegistry.instance().atlas_numel();
            
            % reshape
            GLMmask  = reshape(GLMmask,  [1, N3D]); % 1 x 147456 single 
            featdata = reshape(featdata, [1, N3D]); 
            funcdata = reshape(funcdata, [1, N3D]);             
            
            % select GLMmask
            featdata_glm = featdata(logical(GLMmask)); % [0,1]
            funcdata_glm = funcdata(logical(GLMmask)); % 1 x 65549
            funcdata_glm = (1 - funcdata_glm)/2;       % [0 1] measure of anomaly
            
            % prepare for roc()
            featdata_glm = [featdata_glm; 1 - featdata_glm];
            funcdata_glm = [funcdata_glm; 1 - funcdata_glm];
            [tpr,fpr,thresh] = roc(featdata_glm, funcdata_glm);
            %figure; plot(tpr{1}); title('TPR')
            %figure; plot(fpr{1}); title('FPR')
            %figure; plot(thresh{1}); title('thresholds')
            figure; plotroc(featdata_glm, funcdata_glm)            
        end
        function vec = removeInfratentorial(vec)
            if ndims(vec) < 3 %#ok<ISMAT>
                vec = reshape(vec, [48 64 48]);
            end
            vec(:,:,1:16) = 0;
            if ~isempty(getenv('DEBUG'))
                ic2 = mlfourd.ImagingContext2(vec); ic2.fsleyes
            end
        end
    end

	methods 
        
        %% GET
        
        function g = get.atlas_1d(~)
            g = mlperceptron.PerceptronRegistry.read_MNI152_T1_0p5mm();
        end
        
        %%
		
        function afc_corr = hotspot_corr(this, varargin)
            %% HOTSPOT_CORR
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
            addParameter(ip, 'feature_filename', 'afc.4dfp.hdr', @isfile)
            addParameter(ip, 'load_afc', true, @islogical)
            addParameter(ip, 'Nsigma', 2, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.patientdir = ipr.patientdir;
            this.patientid = ipr.patientid;                        
            workdir = fileparts(ipr.afc_filename);
            pwd0 = pushd(workdir);
            
            % retrieve afc \in R^3; threshold to mask
            afc = mlfourd.ImagingContext2(ipr.feature_filename);
            mu = afc.dipmean();
            sigma = afc.dipstd();
            afc_uthr = afc.uthresh(mu - ipr.Nsigma*sigma);            
            
            % retrieve bold \in R^{3+1}; afc_dyn := bold in R using afc mask
            bold = mlperceptron.PerceptronRelease.bold_frames_mat_to_ImagingContext(ipr.patientdir);            
            bold_img = bold.fourdfp.img;
            Nxyz = [size(bold_img, 1) size(bold_img, 2) size(bold_img, 3)];
            Nt = size(bold_img, 4);
            afc_bool = afc_uthr.logical;
            afc_dyn = zeros(Nt, 1);
            for t = 1:Nt
                bold_frame = bold_img(:,:,:,t);
                afc_dyn(t) = mean(bold_frame(afc_bool));
            end 
            
            % correlate afc_dyn in R with bold in R^{1+3}
            bold_img = reshape(bold_img, [prod(Nxyz) Nt]);
            c = corr(afc_dyn, bold_img');
            c = reshape(c, Nxyz);
            afc_corr = mlfourd.ImagingContext2(c, 'mmppix', [3 3 3], 'filename', 'afc_corr.4dfp.hdr');            
            
            popd(pwd0)
        end   
        function [sim,featifc,funcifc] = makeDice(this, varargin)
            %% MAKEDICE
            %  Preconditions:
            %  -- completion of SL_fMRI_initialization
            %  -- completion of mlperceptron.PerceptronRelease factory
            %  @param patientdir is a folder.
            %  @param patientid is char.
            %  @param afc_filename is char.
            %  @param Nframes is numeric.
            %  @param feature is a filename, NIfTI or 4dfp.
            %  @returns similarity
            %  @returns ImagingFormatContext
            
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
            this.feature = this.imgread(ipr.feature_filename);
            fthresh = dipmax(this.feature)/2;
            this.feature = single(this.feature > fthresh);
            [sim,featifc,funcifc] = this.calcdice(this.removeInfratentorial(this.feature), ...
                                                  this.removeInfratentorial(this.afc_map));
            %saveFigures(workdir)
            popd(pwd0)
        end
        function this = makeDifferenceMap(this, varargin)
            %% MAKEDIFFERENCEMAP
            %  Preconditions:
            %  -- completion of SL_fMRI_initialization
            %  -- completion of mlperceptron.PerceptronRelease factory
            %  @param patientdir is a folder.
            %  @param patientid is char.
            %  @param afc_filename is char.
            %  @param Nframes is numeric.
            %  @returns instance of mlafc.AFC.  
            %  @returns inputParser.Results.       
            
            [this,ipr] = this.makeSearchlightMap(varargin{:});
            
            workdir = fullfile(ipr.patientdir, myfileprefix(ipr.afc_filename));
            mkdir(workdir)
            pwd0 = pushd(workdir);
            save(ipr.afc_filename, 'this')            
            figure
            this.mapOverlay(this.atlas_1d', this.afc_map)              
            saveFigures(workdir)
            popd(pwd0)
        end     
        function this = makeResectionMap(this, varargin)
            %% MAKERESECTIONMAP
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
            addParameter(ip, 'load_afc', false, @islogical)
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
            this.feature = this.imgread(ipr.feature_filename);
            this.feature(this.feature > 0) = 1;
            this.mapOverlay2(this.read_T1_333', this.afc_map, this.feature)                
            saveFigures(workdir)
            popd(pwd0)
        end
        function this = makeROC(this, varargin)
            %% MAKEROC
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
            this.feature = this.imgread(ipr.feature_filename);
            fthresh = dipmax(this.feature)/2;
            this.feature = single(this.feature > fthresh);
            this.plotroc(this.removeInfratentorial(this.feature), ...
                         this.removeInfratentorial(this.afc_map))
            saveFigures(workdir)
            popd(pwd0)
        end    
        function [this,ipr] = makeSoftmax(this, varargin)
            %% MAKESOFTMAX of dissimilarity
            %  Preconditions:
            %  -- completion of SL_fMRI_initialization
            %  -- completion of mlperceptron.PerceptronRelease factory
            %  @param patientdir is a folder.
            %  @param patientid is char.
            %  @param afc_filename is char.
            %  @param Nframes is numeric.
            %  @returns instance of mlafc.AFC.
            %  @returns inputParser.Results.
            %  @returns softmax in product cast as ImagingContext2.
            
            if this.registry.tanh_sandwich
                [this,ipr] = this.makeSoftmax_tanh_atanh__(varargin{:});
            else
                [this,ipr] = this.makeSoftmax__(varargin{:});
            end
        end
        function [this,ipr] = makeSoftmax_tanh_atanh__(this, varargin)
            import mlperceptron.PerceptronFromFormat
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
            
            %% using results of mlperceptron.PerceptronRelease factory                   
            
            PerceptronFromFormat.createProbMaps(this.patientdir, this.patientid)
            
            %% 

            try
                load(fullfile(this.patientdir, 'Perceptron', sprintf('%s%s', this.patientid, this.registry.perceptron_resid_mat)), 'bold_frames')
            catch ME
                handwarning(ME)
                load(fullfile(this.patientdir, 'Perceptron', sprintf('%s%s', this.patientid, this.registry.perceptron_uout_resid_mat)), 'bold_frames')
            end
            % bold_frames is 2307 x 65549 single, time-samples x parenchyma voxels; N_times is variable across studies
            
            if ~isempty(ipr.Nframes)
                bold_frames = bold_frames(1:ipr.Nframes, :);
            end 
            fprintf('makeSearchlightMap.bold_frames.size -> %s\n', mat2str(size(bold_frames)))
            
            similarity__ = @this.similarity;
            sv__ = this.sv;
            Nsv__ = length(sv__);
            GMmsk_for_glm__ = this.GMmsk_for_glm_;
            found_GMmsk_for_glm = find(this.GMmsk_for_glm_); % 18611 x 1
            N_GMmask_for_glm = length(this.GMmsk_for_glm_);
            
            sl_fc = zeros(Nsv__, sum(GMmsk_for_glm__), class(bold_frames));
            parfor isv = 1:Nsv__ % EXPENSIVE without prealloc; 1:2825
                sphere_vox = sv__{isv}; % 13 x 1
                bf_sv = nanmean(bold_frames(:, sphere_vox), 2); %#ok<PFBNS>
                sl_fc(isv, :) = similarity__(bf_sv, bold_frames(:, found_GMmsk_for_glm), 'kind', 'ch');
                % sl_fc is 2825 x 18611, |spheres| x |found gray-matter voxels|
            end
            
            r = this.similarity(sl_fc, this.sl_fc_mean_); % 1 x 2825 <= 2825 x 18611, 2825 x 18611            
            r_glm = zeros(Nsv__, N_GMmask_for_glm); % 1 x 65549
            count = zeros(Nsv__, N_GMmask_for_glm); % 1 x 65549            
            parfor isv = 1:Nsv__               
                sphere_vox = sv__{isv};
                r_glm_inc = zeros(1, N_GMmask_for_glm); 
                r_glm_inc(sphere_vox) = atanh(r(isv)); 
                r_glm(isv,:) = r_glm_inc; 
                count_inc = zeros(1, N_GMmask_for_glm);
                count_inc(sphere_vox) = 1;
                count(isv,:) = count_inc;               
            end            
            energy_pi = tanh(sum(r_glm,1)./sum(count,1)); % 1 x 65549
            prob_pi = exp(-energy_pi);
            
            probs = exp(-this.sl_fc_gsp_); % 100 x 65549
            probs(size(probs,1)+1,:) = prob_pi;
            this.afc_map = prob_pi ./ sum(probs, 1); 
            this.afc_map = this.maskedVecToFullVec(this.afc_map); % 1 x 147456
        end
        function [this,ipr] = makeSoftmax__(this, varargin)
            import mlperceptron.PerceptronFromFormat
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
            
            %% using results of mlperceptron.PerceptronRelease factory                   
            
            PerceptronFromFormat.createProbMaps(this.patientdir, this.patientid)
            
            %% 

            try
                load(fullfile(this.patientdir, 'Perceptron', sprintf('%s%s', this.patientid, this.registry.perceptron_resid_mat)), 'bold_frames')
            catch ME
                handwarning(ME)
                load(fullfile(this.patientdir, 'Perceptron', sprintf('%s%s', this.patientid, this.registry.perceptron_uout_resid_mat)), 'bold_frames')
            end
            % bold_frames is 2307 x 65549 single, time-samples x parenchyma voxels; N_times is variable across studies
            
            if ~isempty(ipr.Nframes)
                bold_frames = bold_frames(1:ipr.Nframes, :);
            end 
            fprintf('makeSearchlightMap.bold_frames.size -> %s\n', mat2str(size(bold_frames)))
            
            similarity__ = @this.similarity;
            sv__ = this.sv;
            Nsv__ = length(sv__);
            GMmsk_for_glm__ = this.GMmsk_for_glm_;
            found_GMmsk_for_glm = find(this.GMmsk_for_glm_); % 18611 x 1
            N_GMmask_for_glm = length(this.GMmsk_for_glm_);
            
            sl_fc = zeros(Nsv__, sum(GMmsk_for_glm__), class(bold_frames));
            parfor isv = 1:Nsv__ % EXPENSIVE without prealloc; 1:2825
                sphere_vox = sv__{isv}; % 13 x 1
                bf_sv = nanmean(bold_frames(:, sphere_vox), 2); %#ok<PFBNS>
                sl_fc(isv, :) = similarity__(bf_sv, bold_frames(:, found_GMmsk_for_glm), 'kind', 'ch');
                % sl_fc is 2825 x 18611, |spheres| x |found gray-matter voxels|
            end
            
            r = this.similarity(sl_fc, this.sl_fc_mean_); % 1 x 2825 <= 2825 x 18611, 2825 x 18611            
            r_glm = zeros(Nsv__, N_GMmask_for_glm); % 1 x 65549
            count = zeros(Nsv__, N_GMmask_for_glm); % 1 x 65549            
            parfor isv = 1:Nsv__               
                sphere_vox = sv__{isv};
                r_glm_inc = zeros(1, N_GMmask_for_glm); 
                r_glm_inc(sphere_vox) = r(isv); 
                r_glm(isv,:) = r_glm_inc; 
                count_inc = zeros(1, N_GMmask_for_glm);
                count_inc(sphere_vox) = 1;
                count(isv,:) = count_inc;               
            end            
            energy_pi = sum(r_glm,1)./sum(count,1); % 1 x 65549
            prob_pi = exp(-energy_pi);
            
            probs = exp(-this.sl_fc_gsp_); % 100 x 65549
            probs(size(probs,1)+1,:) = prob_pi;
            this.afc_map = prob_pi ./ sum(probs, 1); 
            this.afc_map = this.maskedVecToFullVec(this.afc_map); % 1 x 147456
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
            
            import mlperceptron.PerceptronFromFormat
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
            
            %% using results of mlperceptron.PerceptronRelease factory                   
            
            PerceptronFromFormat.createProbMaps(this.patientdir, this.patientid)
            
            %% 

            try
                load(fullfile(this.patientdir, 'Perceptron', sprintf('%s%s', this.patientid, this.registry.perceptron_resid_mat)), 'bold_frames')
            catch ME
                handwarning(ME)
                load(fullfile(this.patientdir, 'Perceptron', sprintf('%s%s', this.patientid, this.registry.perceptron_uout_resid_mat)), 'bold_frames')
            end
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
                % sl_fc is 2825 x 18611, |spheres| x |found gray-matter voxels|
                
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
        function plotIntersubjectVariability(this)
            
            import mlpark.SearchLight
            
            addpath(genpath('Z:\shimony\Park\Matlab'))
            tgtatlas = 'Z:\shimony\Park\atlas\Targets\MNI152_711-2B_333.4dfp.img';
            tgtanatData = this.imgread(tgtatlas);
            
            %% MLP AFC
            numsub = 150;
            MLP_GTM_150 = this.gtm_y; % GTM_y(:,:,1:numsub);
            
            ind = 1;
            for sub1 = 1:numsub
                for sub2 = sub1+1:numsub
                    
                    MLP_GTM_sub1 = MLP_GTM_150(:,:,sub1);
                    MLP_GTM_sub2 = MLP_GTM_150(:,:,sub2);
                    
                    % rmse
                    sq = (MLP_GTM_sub1 - MLP_GTM_sub2).^2;
                    m = sum(sq')/8; %#ok<UDIM>
                    r = sqrt(m);
                    mlp_rmse_con150(:,ind) = r; 
                    
                    % cc
                    mlp_cc(:,ind) = SearchLight.corr_v2(MLP_GTM_sub1, MLP_GTM_sub2); 
                    
                    ind = ind+1;
                end
            end
            
            %%
            E_mean = mean(mlp_rmse_con150'); %#ok<UDIM>
            
            intersubject_var = zeros(1,147456);
            intersubject_var(this.glmmsk_indices_) = E_mean;
            
            figure; 
            SearchLight.MapOverlay(tgtanatData, intersubject_var,'caxvals',[0 0.2],'ttlvec', "Intersubject Variability: GSP 100")
            
            %%
            cc_mean = mean(mlp_cc'); %#ok<UDIM>
            
            intersubject_var = zeros(1,147456);
            intersubject_var(this.glmmsk_indices_) = cc_mean;
            
            figure; 
            SearchLight.MapOverlay(tgtanatData, intersubject_var,'caxvals',[0 1],'ttlvec', "Intersubject Variability: GSP 100")
        end
        function p = product(this)
            perreg = mlperceptron.PerceptronRegistry.instance();
            [d1,d2,d3] = perreg.atlas_dims;
            img = reshape(this.afc_map, [d1 d2 d3]);
            img(isnan(img)) = 0;
            warning('mlafc:RuntimeWarning', ...
                'AFC.product() is flipping axes 1 & 2 to correct flips from mlperceptron.Fourdfp.Read4dfp')
            img = flip(flip(img, 1), 2);
            
            ifc = mlfourd.ImagingFormatContext( ...
                fullfile(getenv('REFDIR'), '711-2B_333.4dfp.hdr'));
            ifc.img = img;
            ifc.filepath = pwd;
            ifc.fileprefix = 'AFC_product';
            p = mlfourd.ImagingContext2(ifc);
        end
        function [img,frames,voxelsize] = read_T1_333(this)
            import mlperceptron.*
            [img,frames,voxelsize] = Fourdfp.read_4dfpimg_v1( ...
                fullfile(this.patientdir, 'atlas', [this.patientid '_mpr_n1_333_t88.4dfp.img']));
        end  
        function cut_corr = resection_corr(this, varargin)
            %% RESECTION_CORR
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
            addParameter(ip, 'Nsigma', 2, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.patientdir = ipr.patientdir;
            this.patientid = ipr.patientid;                        
            workdir = fileparts(ipr.afc_filename);
            pwd0 = pushd(workdir);
            
            % retrieve afc \in R^3; threshold to mask
            cut = mlfourd.ImagingContext2(sprintf('../%s_Segmentation.4dfp.hdr', this.patientid));
            mu = cut.dipmean();
            cut_thr = cut.thresh(mu);
            
            % retrieve bold \in R^{3+1}; cut_dyn := bold in R using afc mask
            bold = mlperceptron.PerceptronRelease.bold_frames_mat_to_ImagingContext(ipr.patientdir);            
            bold_img = bold.fourdfp.img;
            Nxyz = [size(bold_img, 1) size(bold_img, 2) size(bold_img, 3)];
            Nt = size(bold_img, 4);
            cut_bool = cut_thr.logical;
            cut_dyn = zeros(Nt, 1);
            for t = 1:Nt
                bold_frame = bold_img(:,:,:,t);
                cut_dyn(t) = mean(bold_frame(cut_bool));
            end 
            
            % correlate cut_dyn in R with bold in R^{1+3}
            bold_img = reshape(bold_img, [prod(Nxyz) Nt]);
            c = corr(cut_dyn, bold_img');
            c = reshape(c, Nxyz);
            cut_corr = mlfourd.ImagingContext2(c, 'mmppix', [3 3 3], 'filename', 'resection_corr.4dfp.hdr');            
            
            popd(pwd0)
        end      
        function s = similarity(this, x, y, varargin)
            ip = inputParser;
            addParameter(ip, 'kind', this.similarityKind, @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            switch lower(ipr.kind)
                case 'ch'
                    s = mlperceptron.PerceptronRelease.corr_new(x, y);
                case 'kp'
                    s = mlpark.SearchLight.corr_kp(x, y);
                case 'kl'
                    s = 1 - tanh(mlafc.AFC.kldiv(x, y)); % dissimilarity -> similarity
                case 'js'
                    s = 1 - tanh(mlafc.AFC.kldiv(x, y, 'js')); % dissimilarity -> similarity
                otherwise
                    error('mlafc:RuntimeError', 'AFC.similarity.ipr.kind->%s', ipr.kind)
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
            addParameter(ip, 'make_sl_fc_mean', true, @islogical)
            addParameter(ip, 'make_sl_fc_gsp', true, @islogical)
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
%                    sl_fcz = atanh(sl_fc);
                    sl_fcz = sl_fc;
                    sl_fcz_sum = sl_fcz_sum + sl_fcz;
                    
                    fprintf('%s processed: %f\n', num2str(refnum), toc)
                    clear bold_frames sl_fc sl_fcz                    
                end
                
                % average
                sl_fcz_mean = sl_fcz_sum / this.registry_.ref_count;
                % tanh <-> inverse Fisher z-transform
%                sl_fc_mean = tanh(sl_fcz_mean);
                sl_fc_mean = sl_fcz_mean;
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
%                        temp(this.sv_{ind}) = atanh(r(ind)); 
                        temp(this.sv_{ind}) = r(ind); 
                        counttemp(this.sv_{ind}) = 1;
                        r_glm = r_glm + temp; 
                        count = count + counttemp;                        
                    end
                    
                    % tanh <-> inverse Fisher z-transform
%                    sl_fc_gsp(refnum,:) = tanh(r_glm./count);
                    sl_fc_gsp(refnum,:) = r_glm./count;
                    fprintf('%s processed: %f\n', num2str(refnum), toc)                    
                end
                save(this.registry_.sl_fc_gsp_mat, 'sl_fc_gsp')
                this.sl_fc_gsp_ = sl_fc_gsp;
                clear('sl_fc_gsp')
                popd(pwd0)
            end
        end
        function this = SL_initialization(this, varargin)
            %% SL_INITIALIZATION
            %  @return this.sl_fc_mean_: 
            %          from %% SL-fMRI mean map using 100 GSP subjects
            %          2825 (# spheres) x 18611 (# grey matter voxels)
            %  @return this.sl_fc_gsp_:
            %          from %% SL-fMRI correlation maps of 100 GSP subjects
            %          100 (# subjects) x 65549 (# whole brain voxels)
            %          this variable represents inter-subject variability across subjects in each voxel  
            
            % 04/18/19 KP
            % 06/06/19 JJL          
            %
            % If you want to change the sphere radius or grid spacing, you should run
            % this code again to create a new SL-fMRI "ground truth map"
            
            if this.registry.tanh_sandwich
                this = this.SL_tanh_atanh_initialization__(varargin{:});
            else
                this = this.SL_initialization__(varargin{:});
            end
        end           
        function this = SL_tanh_atanh_initialization__(this, varargin)
            import mlperceptron.*
            import mlperceptron.PerceptronRelease;
            import mlpark.SearchLight;
            
            reg = this.registry;
            pwd0 = pushd(reg.gtm500_dir); 
            
            gtm500_dir = reg.gtm500_dir;
            ref_resid_mat = reg.ref_resid_mat;
            similarity__ = @this.similarity;            
            found_GMmsk_for_glm = find(this.GMmsk_for_glm_); % 18611 x 1
            N_GMmsk_for_glm = length(this.GMmsk_for_glm_);
            sv__ = this.sv;
            Nsv__ = length(sv__);            
            
            %% make SL mean map using 100 GSP subjects
          
            sl_fc_cache = zeros(length(this.sv), length(found_GMmsk_for_glm)); 

            for refnum = 1:reg.ref_count
                refid = reg.gtm500_ids{refnum};
                ld1 = load(fullfile(gtm500_dir, ...
                           refid, ...
                           'Perceptron', ...
                           sprintf('%s%s', refid, ref_resid_mat)),'bold_frames'); % bold_frames has 245 x 65549  
                bold_frames = ld1.bold_frames;
                           
                parfor isv = 1:Nsv__                     
                    sphere_vox = sv__{isv}; % 13 x 1
                    bf_sv = mean(bold_frames(:, sphere_vox),2); %#ok<PFBNS> % 245 x 1
                    sl_fc_cache(isv,:) = sl_fc_cache(isv,:) + ...
                        atanh(similarity__(bf_sv, bold_frames(:,found_GMmsk_for_glm), 'kind', 'ch')); % SL_FC
                        % 2825 x 18611
                end                    
            end
            sl_fc_mean = tanh(sl_fc_cache./reg.ref_count);
            save(reg.sl_fc_mean_mat, 'sl_fc_mean')
            
            %% make SL correlation maps of 100 GSP subjects
                   
            sl_fc_cache = zeros(length(this.sv), length(found_GMmsk_for_glm));              
            sl_fc_gsp = zeros(reg.ref_count, N_GMmsk_for_glm); 
            
            for refnum = 1:reg.ref_count
                refid = reg.gtm500_ids{refnum};
                ld1 = load(fullfile(gtm500_dir, ...
                           refid, ...
                           'Perceptron', ...
                           sprintf('%s%s', refid, ref_resid_mat)),'bold_frames'); % bold_frames has 245 x 65549  
                bold_frames = ld1.bold_frames;

                parfor isv = 1:Nsv__                     
                    sphere_vox = sv__{isv}; % 13 x 1
                    bf_sv = mean(bold_frames(:, sphere_vox),2); %#ok<PFBNS> % 245 x 1
                    sl_fc_cache(isv,:) = ...
                        atanh(similarity__(bf_sv, bold_frames(:,found_GMmsk_for_glm), 'kind', 'ch')); % SL_FC
                        % 2825 x 18611
                end

                r = similarity__(sl_fc_cache, sl_fc_mean); % contract gm voxels, return 1 x length(sv__)                                                    
                r_glm = zeros(Nsv__, N_GMmsk_for_glm);
                count = zeros(Nsv__, N_GMmsk_for_glm);
                parfor isv = 1:Nsv__ % spheres will overlap for stride < sphere diameter
                    sphere_vox = sv__{isv};
                    r_glm_inc = zeros(1, N_GMmsk_for_glm);
                    r_glm_inc(sphere_vox) = atanh(r(isv));
                    r_glm(isv,:) = r_glm_inc;
                    count_inc = zeros(1, N_GMmsk_for_glm);
                    count_inc(sphere_vox) = 1;
                    count(isv,:) = count_inc;
                end % f_{sigma -> xi^pi}
                sl_fc_gsp(refnum,:) = tanh(sum(r_glm,1)./sum(count,1));
            end
            save(reg.sl_fc_gsp_mat, 'sl_fc_gsp')
            
            %% finalize          
            
            this.sl_fc_mean_ = sl_fc_mean;
            this.sl_fc_gsp_ = sl_fc_gsp;
            clear('sl_fc_mean')
            clear('sl_fc_gsp')
            
            popd(pwd0)
        end 
        function this = SL_initialization__(this, varargin)
            import mlperceptron.*
            import mlperceptron.PerceptronRelease;
            import mlpark.SearchLight;
            
            reg = this.registry;
            pwd0 = pushd(reg.gtm500_dir); 
            
            gtm500_dir = reg.gtm500_dir;
            ref_resid_mat = reg.ref_resid_mat;
            similarity__ = @this.similarity;            
            found_GMmsk_for_glm = find(this.GMmsk_for_glm_); % 18611 x 1
            N_GMmsk_for_glm = length(this.GMmsk_for_glm_);
            sv__ = this.sv;
            Nsv__ = length(sv__);            
            
            %% make SL mean map using 100 GSP subjects
          
            sl_fc_cache = zeros(length(this.sv), length(found_GMmsk_for_glm)); 

            for refnum = 1:reg.ref_count
                refid = reg.gtm500_ids{refnum};
                ld1 = load(fullfile(gtm500_dir, ...
                           refid, ...
                           'Perceptron', ...
                           sprintf('%s%s', refid, ref_resid_mat)),'bold_frames'); % bold_frames has 245 x 65549  
                bold_frames = ld1.bold_frames;
                           
                parfor isv = 1:Nsv__                     
                    sphere_vox = sv__{isv}; % 13 x 1
                    bf_sv = mean(bold_frames(:, sphere_vox),2); %#ok<PFBNS> % 245 x 1
                    sl_fc_cache(isv,:) = sl_fc_cache(isv,:) + ...
                        similarity__(bf_sv, bold_frames(:,found_GMmsk_for_glm), 'kind', 'ch'); % SL_FC
                        % 2825 x 18611
                end                    
            end
            sl_fc_mean = sl_fc_cache./reg.ref_count;
            save(reg.sl_fc_mean_mat, 'sl_fc_mean')
            
            %% make SL correlation maps of 100 GSP subjects
                   
            sl_fc_cache = zeros(length(this.sv), length(found_GMmsk_for_glm));              
            sl_fc_gsp = zeros(reg.ref_count, N_GMmsk_for_glm); 
            
            for refnum = 1:reg.ref_count
                refid = reg.gtm500_ids{refnum};
                ld1 = load(fullfile(gtm500_dir, ...
                           refid, ...
                           'Perceptron', ...
                           sprintf('%s%s', refid, ref_resid_mat)),'bold_frames'); % bold_frames has 245 x 65549  
                bold_frames = ld1.bold_frames;

                parfor isv = 1:Nsv__                     
                    sphere_vox = sv__{isv}; % 13 x 1
                    bf_sv = mean(bold_frames(:, sphere_vox),2); %#ok<PFBNS> % 245 x 1
                    sl_fc_cache(isv,:) = ...
                        similarity__(bf_sv, bold_frames(:,found_GMmsk_for_glm), 'kind', 'ch'); % SL_FC
                        % 2825 x 18611
                end

                r = similarity__(sl_fc_cache, sl_fc_mean); % contract gm voxels, return 1 x length(sv__)                                                    
                r_glm = zeros(Nsv__, N_GMmsk_for_glm);
                count = zeros(Nsv__, N_GMmsk_for_glm);
                parfor isv = 1:Nsv__ % spheres will overlap for stride < sphere diameter
                    sphere_vox = sv__{isv};
                    r_glm_inc = zeros(1, N_GMmsk_for_glm);
                    r_glm_inc(sphere_vox) = r(isv);
                    r_glm(isv,:) = r_glm_inc;
                    count_inc = zeros(1, N_GMmsk_for_glm);
                    count_inc(sphere_vox) = 1;
                    count(isv,:) = count_inc;
                end % f_{sigma -> xi^pi}
                sl_fc_gsp(refnum,:) = sum(r_glm,1)./sum(count,1);
            end
            save(reg.sl_fc_gsp_mat, 'sl_fc_gsp')
            
            %% finalize          
            
            this.sl_fc_mean_ = sl_fc_mean;
            this.sl_fc_gsp_ = sl_fc_gsp;
            clear('sl_fc_mean')
            clear('sl_fc_gsp')
            
            popd(pwd0)
        end
        
 		function this = KaysAFC(varargin)
 			%% KAYSAFC
 			%  @param .

 			this = this@AFC(varargin{:});
            
            this = sv_initialization(this);
            this = load(this);
 		end
    end 
    
    %% PROTECTED    
    
    methods (Static, Access = protected)
        function im  = snapImage(image, snapsize, snapind, plane) 
            %% Supports only 333

            [Nx,Ny] = mlafc.AFCRegistry.instance().atlas_dims();
            switch plane
                case 'axial'
                    Nx = 48; Ny = 64;
                    snapdefind = reshape(1:48, [12 4])';
                case 'coronal'
                    Nx = 48; Ny = 48;
                    snapdefind = reshape(1:64, [8 8])';
            end

            xsi = 1;
            for si = 1:length(snapind)                    
                if si == 1
                    xsi = 1;
                else
                    if xsi < snapsize(2)
                        xsi = xsi+1;
                    else
                        xsi = 1;
                    end
                end
                ysi = ceil(si/snapsize(2));

                [i, j] = find(snapdefind == snapind(si));
                imtemp = image(Ny*(i-1)+1:Ny*i, Nx*(j-1)+1:Nx*j);
                im(Ny*(ysi-1)+1:Ny*ysi, Nx*(xsi-1)+1:Nx*xsi) = imtemp;                    
            end
        end
    end
    
    methods (Access = protected)  
        function map = diff_map(this, sl_fmri_pat)
            %% DIFF_MAP
            %  @param sl_fmri_pat is vec, numel(vec) == 65549.
            %  @returns map of difference from this.sl_fc_gsp_ that is 1 x 147456.
            
            sl_fc_gsp_mean = nanmean(this.sl_fc_gsp_); % 1 x 65549            
            afc_glmmap = sl_fmri_pat - sl_fc_gsp_mean; % 1 x 65549  
            map = this.maskedVecToFullVec(afc_glmmap);
        end    
        function g = gtm_y(this)               
            ld = load(fullfile(this.registry.gtm500_dir, 'GTM_y.mat'), 'GTM_y');
            g = ld.GTM_y;
        end
        function this = load(this)
            if isfile(this.registry.sl_fc_mean_mat)
                mat = load(this.registry.sl_fc_mean_mat);
                this.sl_fc_mean_ = mat.sl_fc_mean;
            end
            if isfile(this.registry.sl_fc_gsp_mat)
                mat = load(this.registry.sl_fc_gsp_mat);
                this.sl_fc_gsp_ = mat.sl_fc_gsp;
            end
        end  
        function this = sv_initialization(this)
            %% SV_INITIALIZATION
            %  @return this.sv_.
            
            % 04/18/19 KP
            % 06/06/19 JJL   
            
            %% Store sphere voxels indices
            
            this.sv_ = {};
            [Nx,Ny,Nz] = this.registry.atlas_dims;
            glmmsk_3d = reshape(this.glmatl_, [Nx, Ny, Nz]);
            R = this.registry.sphere_radius;
            I = this.registry.grid_spacing;
            isv = 1;
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
                        if length(sphere_vox) >= this.registry.min_num_vox
                            this.sv_{isv} = sphere_vox;
                            isv = isv+1;
                        end
                    end
                end
            end
        end
        function vec = maskedVecToFullVec(this, vec0)
            %% MASKEDVECTOFULLVEC
            %  @param   vec0 is Nt x 65549 
            %  @returns vec  is Nt x 147456
            
            assert(isnumeric(vec0) && length(vec0) == 65549)
            
            vec = zeros(1, this.registry.atlas_numel()); 
            vec(this.glmmsk_indices_) = vec0; 
        end
        function map = threshed_afc_map(this, varargin)
            %% THRESHED_AFC_MAP
            %  @param sl_fmri_pat is vec, numel(vec) <= 65549.
            %  @param sd_thr is z-threshold.
            %  @returns map, a vec that is 1 x 147456.
            
            ip = inputParser;
            addRequired(ip, 'sl_fmri_pat', @isnumeric)
            addOptional(ip, 'sd_thr', 3, @isnumeric);
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            sl_fc_gsp_mean = nanmean(this.sl_fc_gsp_);              % 1 x 65549
            sl_fc_gsp_std  = nanstd(this.sl_fc_gsp_);               % 1 x 65549
            sl_fc_gsp_thr  = sl_fc_gsp_mean - ipr.sd_thr*sl_fc_gsp_std; % 1 x 65549

            afc_vox             = find(ipr.sl_fmri_pat < sl_fc_gsp_thr);     % 1 x 56281
            afc_glmmap          = zeros(1, length(this.GMmsk_for_glm_)); % 1 x 65549
            afc_glmmap(afc_vox) = 1;                                     % 1 x 65549            
            map = this.maskedVecToFullVec(afc_glmmap);
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

