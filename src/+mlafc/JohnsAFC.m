classdef JohnsAFC < mlafc.AFC
	%% JOHNSAFC  

	%  $Revision$
 	%  was created 04-Apr-2021 12:25:45 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlafc/src/+mlafc.
 	%% It was developed on Matlab 9.9.0.1592791 (R2020b) Update 5 for MACI64.  Copyright 2021 John Joowon Lee.
 	
    properties (Constant)
        N_BOLD = 65549
    end
    
	properties (Dependent)
 		GMctx % ImagingContext2
        GMToYeo7_xfm
        BigBrain300
        Yeo7
        Yeo17
    end
    
    methods (Static)
        function acc = acorrcoef(A, B)
            %% ACORRCOEF builds the corrcoef (Ns x Nv) for A (Nt x Ns) and B (Nt x Nv).
            %  Ns ~ num spheres.
            %  Nv ~ num voxels.
            
            [Nt,Ns] = size(A);
            [Nt_,Nv] = size(B);
            assert(Nt == Nt_)
            
            muA = mean(A, 1); % 1 x Ns
            muB = mean(B, 1); % 1 x Nv
            
            dA = A - muA(ones(Nt,1),:); % Nt x Ns
            dB = B - muB(ones(Nt,1),:); % Nt x Nv
            
            cov = dA'*dB; % Ns x Nv
            
            sigmaA = sqrt(sum(dA.^2, 1))'; % Ns x 1
            sigmaB = sqrt(sum(dB.^2, 1));  % 1  x Nv
            
            acc = cov ./ (sigmaA(:,ones(1,Nv)) .* sigmaB(ones(Ns,1),:)); % Ns x Nv
        end
        function arr = glmmskArrToFullArr(arr0)
            %% GLMMSKARRTOFULLARR
            %  @param required arr0 is Nt x 65549.
            %  @return arr is Nt x (48*64*48).
            
            assert(isnumeric(arr0))
            assert(size(arr0,2) == 65549)
            assert(ismatrix(arr0))            
            
            glmatl_ = mlperceptron.PerceptronRegistry.read_glm_atlas_mask(); % 147456 x 1
            glmatl_(find(glmatl_)) = 1;   %#ok<FNDSB>
            found_glmatl_ = find(glmatl_); % 65549 x 1
            
            Nt = size(arr0,1);
            arr = zeros(Nt, mlafc.AFCRegistry.instance().atlas_numel()); % Nt x 147456
            if 1 == Nt
                arr(found_glmatl_') = arr0;
                return
            end            
            for t = 1:Nt
                arr(t,found_glmatl_') = arr0(t,:);
            end
        end     
        function arr = GMmskArrToFullArr(arr0)
            %% GMMSKARRTOFULLARR
            %  @param required arr0 is Npheres x 18611.
            %  @return arr is Nspheres x (48*64*48).
            
            import mlperceptron.PerceptronRegistry
            
            assert(isnumeric(arr0))
            assert(size(arr0,2) == 18611)
            assert(ismatrix(arr0)) 
            
%             glmatl_ = PerceptronRegistry.read_glm_atlas_mask(); % 147456 x 1
%             glmatl_(find(glmatl_)) = 1; %#ok<FNDSB>
%             GMctx_ = PerceptronRegistry.read_N21_aparc_aseg_GMctx(); % 147456 x 1
%             GMmsk_for_glm_ = GMctx_(find(glmatl_)); %#ok<FNDSB>
%             GMmsk_indices_ = find(GMmsk_for_glm_); % 18611 x 1            
            
            GMctx_ = PerceptronRegistry.read_N21_aparc_aseg_GMctx(); % 147456 x 1
            found_GMctx_ = find(GMctx_); % 18611 x 1  
            
            Nspheres = size(arr0,1);
            arr = zeros(Nspheres, mlafc.AFCRegistry.instance().atlas_numel()); % Nspheres x 147456
            if 1 == Nspheres
                arr(found_GMctx_') = arr0;
                return
            end            
            for t = 1:Nspheres
                arr(t,found_GMctx_') = arr0(t,:);
            end
        end      
        function ic = slfcMatToIC(matname, objname, varargin)
            %% slfcMatToIC
            %  @param matname is a Percptron mat file.
            %  @param objname is the name of the sought object.
            %  @param flip1 is logical.
            %  @param flip2 is logical.
            %  @param refnum is scalar.
            
            ip = inputParser;
            addRequired(ip, 'matname', @isfile)
            addRequired(ip, 'objname', @ischar)
            addParameter(ip, 'flip1', true, @islogical)
            addParameter(ip, 'flip2', true, @islogical)
            addParameter(ip, 'refnum', 1, @isscalar)
            parse(ip, matname, objname, varargin{:})
            ipr = ip.Results;
            
            assert(isfile(matname))
            assert(ischar(objname))
            ld = load(matname, objname);
            obj = ld.(objname);
            if 3 == ndims(obj)
                obj = squeeze(obj(ipr.refnum, :, :));
            end
            obj = mlafc.JohnsAFC.glmmskArrToFullArr(obj);
            obj = obj';
            
            Nspheres = size(obj, 2);
            img = reshape(obj, [48 64 48 Nspheres]);
            img(isnan(img)) = 0;
            
            if ipr.flip1                
                img = flip(img, 1);
            end
            if ipr.flip2
                img = flip(img, 2);
            end
            ifc = mlfourd.ImagingFormatContext( fullfile(getenv('REFDIR'), '711-2B_333.nii.gz'));
            ifc.img = img;
            ifc.filepath = pwd;
            ifc.fileprefix = myfileprefix(matname);
            ic = mlfourd.ImagingContext2(ifc);
        end
        function ic = perceptronMatToIC(matname, objname, varargin)
            %% PERCEPTRONMATTOIC
            %  @param matname is a Percptron mat file.
            %  @param objname is the name of the sought object.
            %  @param flip1 is logical.
            %  @param flip2 is logical.
            
            ip = inputParser;
            addRequired(ip, 'matname', @isfile)
            addRequired(ip, 'objname', @ischar)
            addParameter(ip, 'flip1', true, @islogical)
            addParameter(ip, 'flip2', true, @islogical)
            parse(ip, matname, objname, varargin{:})
            ipr = ip.Results;
            
            assert(isfile(matname))
            assert(ischar(objname))
            ld = load(matname, objname);
            obj = ld.(objname);
            obj = mlafc.JohnsAFC.glmmskArrToFullArr(obj);
            obj = obj';
            
            Nt = size(obj, 2);
            img = reshape(obj, [48 64 48 Nt]);
            img(isnan(img)) = 0;
            
            if ipr.flip1                
                img = flip(img, 1);
            end
            if ipr.flip2
                img = flip(img, 2);
            end
            ifc = mlfourd.ImagingFormatContext( fullfile(getenv('REFDIR'), '711-2B_333.nii.gz'));
            ifc.img = img;
            ifc.filepath = pwd;
            ifc.fileprefix = myfileprefix(matname);
            ic = mlfourd.ImagingContext2(ifc);
        end
    end

	methods        
        
        %% GET
        
        function g = get.GMctx(~)
            g = mlfourd.ImagingContext2( ...
                fullfile(getenv('REFDIR'),  'N21_aparc+aseg_GMctx_on_711-2V_333_avg_zlt0.5_gAAmask_v1_binarized.nii.gz'));
        end
        function g = get.GMToYeo7_xfm(this)
            persistent xfm_
            if isempty(xfm_)
                xfm_ = this.build_GMToYeo7_xfm_();
            end
            g = xfm_;
        end
        function g = get.BigBrain300(~)
            g = mlfourd.ImagingContext2( ...
                fullfile(getenv('REFDIR'), 'BigBrain300', 'BigBrain300_711-2b_allROIs.nii.gz'));
        end
        function g = get.Yeo7(~)
            g = mlfourd.ImagingContext2( ...
                fullfile(getenv('REFDIR'), 'Yeo', 'Yeo2011_7Networks_333.nii.gz'));
        end 
        function g = get.Yeo17(~)
            g = mlfourd.ImagingContext2( ...
                fullfile(getenv('REFDIR'), 'Yeo', 'Yeo2011_17Networks_333.nii.gz'));
        end        
        
        %%
        
        function this = explore_fc(this, varargin)
            %% EXPLORE_FC
            %  cpu time for each refnum ~ 100 s
            %  RAM for each refnum ~ 420 MB
            
            reg = this.registry;
            pwd0 = pushd(reg.gtm500_dir); 
            gtm500_dir = reg.gtm500_dir;
            ref_resid_mat = reg.ref_resid_mat;  
            
            %% make connectivity maps of GSP subjects
                     
            sv__ = this.sv;
            Nsv__ = length(sv__);
            
            sl_fc_accum = zeros(Nsv__, this.N_BOLD); 
            accum = 0;
            for refnum = 1:reg.ref_count
                try
                    refid = reg.gtm500_ids{refnum};
                    ld = load(fullfile(gtm500_dir, ...
                               refid, ...
                               'Perceptron', ...
                               sprintf('%s%s', refid, ref_resid_mat)), 'bold_frames'); % bold_frames has N_t x 65549  
                    bold_frames = ld.bold_frames;

                    sl_fc_gsp_ref = zeros(Nsv__, this.N_BOLD); 
                    acorrcoef = @this.acorrcoef;
                    parfor isv = 1:Nsv__
                        bf_sv = mean(bold_frames(:, sv__{isv}), 2, 'omitnan'); %#ok<PFBNS> % N_t x 1
                        sl_fc_gsp_ref(isv,:) = ...
                            acorrcoef(bf_sv, bold_frames); % N_{sphere_vox} x 65549
                    end
                    save(reg.sl_fc_gsp_ref_mat(refnum), 'sl_fc_gsp_ref')
                    sl_fc_accum = sl_fc_accum + sl_fc_gsp_ref;
                    accum = accum + 1;
                catch ME
                    handwarning(ME)
                end
            end
            
            %% make mean map using GSP subjects
            
            sl_fc_mean = sl_fc_accum/accum;
            
            %% finalize          
            
            save(reg.sl_fc_mean_mat, 'sl_fc_mean')
            this.sl_fc_mean_ = sl_fc_mean;
            clear('sl_fc_mean')
            
            popd(pwd0)
        end  
        function e = energy_fc(this, fc)
            %e = diag(this.acorrcoef(fc, this.sl_fc_mean_))';
            e = mean(this.acorrcoef(fc, this.sl_fc_mean_), 2, 'omitnan')';
        end      
        function [this,ipr] = makeSoftmax(this, varargin)
            %% MAKESOFTMAX of dissimilarity requires completion of explore_fc() which stores 
            %  this.sl_fc_gsp_, this.sl_fc_mean_.
            %  @param patientdir is a folder.
            %  @param patientid is char.
            %  @param Nframes is numeric.  If ~empty(Nframes), inferences uses only specified Nframes.
            %  @returns instance of mlafc.JohnsAFC.
            %  @returns inputParser.Results.
            %  @returns softmax in this.product as ImagingContext2.
            
            import mlpark.*
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'patientdir', @isfolder)
            addRequired( ip, 'patientid', @ischar)
            addParameter(ip, 'Nframes', [], @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;   
            this.patientdir = ipr.patientdir;
            this.patientid = ipr.patientid;
            
            %% per patients, use results of mlperceptron.PerceptronRelease factory                   
            
            mlperceptron.PerceptronFromFormat.createProbMaps(this.patientdir, this.patientid)
            try
                ld = load( ...
                    fullfile(this.patientdir, 'Perceptron', ...
                    sprintf('%s%s', this.patientid, this.registry.perceptron_resid_mat)), ...
                    'bold_frames');
            catch ME
                handwarning(ME)
                ld = load( ...
                    fullfile(this.patientdir, 'Perceptron', ...
                    sprintf('%s%s', this.patientid, this.registry.perceptron_uout_resid_mat)), ...
                    'bold_frames');
            end
            bold_frames = ld.bold_frames; % N_t x N_vxls ~ 2307 x 65549 single
                                          % N_t is variable across studies            
            if ~isempty(ipr.Nframes)
                bold_frames = bold_frames(1:ipr.Nframes, :);                
                fprintf('makeSoftmax.bold_frames.size -> %s\n', mat2str(size(bold_frames)))
            end           
            
            %% per patient, build connectivity map
            
            sv__ = this.sv;
            Nsv__ = length(sv__);
            
            sl_fc = zeros(Nsv__, this.N_BOLD, class(bold_frames)); 
            acorrcoef = @this.acorrcoef;
            parfor isv = 1:Nsv__ % 1:N_{sphere_vox}
                bf_sv = mean(bold_frames(:, sv__{isv}), 2, 'omitnan'); %#ok<PFBNS> % N_t x 1
                sl_fc(isv, :) = ...
                    acorrcoef(bf_sv, bold_frames); % N_{sphere_vox} x 65549
            end
            
            %% project downsampled fibers to base manifold, using mean abs deviation of fiber from mean field of controls
            
            prob = exp(-this.energy_fc(sl_fc)); % 1 x this.N_BOLD

            %% assemble softmax
            
            sum_prob = prob + this.sum_prob_refs(); % init with patient
            this.afc_map = prob ./ sum_prob; % Boltzmann distribution
            this.afc_map = this.glmmskArrToFullArr(this.afc_map); % 1 x 147456
        end
        function sum_prob = sum_prob_refs(this)
            if isfile(this.registry.sl_fc_gsp_sum_prob_mat)
                ld = load(this.registry.sl_fc_gsp_sum_prob_mat);
                sum_prob = ld.sum_prob;
                return
            end
            sum_prob = zeros(1, this.N_BOLD, 'single');
            for refnum = 1:this.registry.ref_count
                try
                    ld = load(this.registry.sl_fc_gsp_ref_mat(refnum));
                    sum_prob = sum_prob + exp(-this.energy_fc(ld.sl_fc_gsp_ref));
                catch ME
                    handwarning(ME)
                end
            end
            save(this.registry.sl_fc_gsp_sum_prob_mat, 'sum_prob')
        end
        function ic = mapOfSpheres(this)
            reg = this.registry;
            glmarr = zeros(1, 65549);
            for isv = 1:length(this.sv)
                glmarr(this.sv{isv}') = isv;
            end
            fullarr = this.glmmskArrToFullArr(glmarr);
            
            img = reshape(fullarr, [48 64 48]);
            img(isnan(img)) = 0;
                          
            img = flip(img, 1);
            img = flip(img, 2);
            ifc = mlfourd.ImagingFormatContext( fullfile(getenv('REFDIR'), '711-2B_333.nii.gz'));
            ifc.img = img;
            ifc.filepath = pwd;
            ifc.fileprefix = sprintf('mapOfSpheres_radius%i_stride%i_N%i', reg.sphere_radius, reg.grid_spacing, reg.ref_count);
            ic = mlfourd.ImagingContext2(ifc);
        end
        
 		function this = JohnsAFC(varargin)
 			%% JOHNSAFC

 			this = this@mlafc.AFC(varargin{:});             
            this.registry_.min_num_vox = 1;
            this.registry_.tanh_sandwich = false;
            this.registry_.tag = '_JohnsAFC';
            
            if isempty(this.sl_fc_mean_)
                this = load(this);
            end
        end        
    end 
    
    %% PROTECTED
    
    properties (Access = protected)
    end
    
    methods (Access = protected)
        function xfm = build_GMToYeo7_xfm_(this)
            
            gm = this.GMctx.nifti;
            gm.img = flip(flip(gm.img, 1), 2);
            gm_vec_ = reshape(gm.img, [1 48*64*48]);
            gm_vec = gm_vec_;
            gm_vec(gm_vec == 0) = nan;
            Ngm = sum(gm_vec_, 'omitnan');
            found_gm = find(gm_vec_);
            
            yeo7 = this.Yeo7.nifti;
            yeo7.img = flip(flip(yeo7.img, 1), 2);
            yeo7_vec = reshape(yeo7.img, [1 48*64*48]);
            yeo7_vec = yeo7_vec .* gm_vec;
            yeo7_vec(yeo7_vec == 0) = 8;
            yeo7_vec(isnan(yeo7_vec)) = 0;
            
            found_yeo = [];
            for rsn = 1:8
                found_yeo = [found_yeo find(yeo7_vec == rsn)]; %#ok<*AGROW>
            end
            assert(numel(found_yeo) == Ngm)
            
            xfm = zeros(Ngm, Ngm);
            for n = 1:Ngm
                xfm(find(found_yeo == found_gm(n)), n) = 1; %#ok<FNDSB>
            end
        end
        function this = load(this)
            if isfile(this.registry.sl_fc_mean_mat)
                mat = load(this.registry.sl_fc_mean_mat);
                this.sl_fc_mean_ = mat.sl_fc_mean;
            end
        end
        function this = sv_initialization(this)
            %% SV_INITIALIZATION
            %  @return this.sv_.
            
            % 04/18/19 KP
            % 06/06/19 JJL   
            
            %% Store sphere voxels indices
            
            glmatl_ = mlperceptron.PerceptronRegistry.read_glm_atlas_mask();
            glmatl_(find(glmatl_)) = 1; %#ok<FNDSB>
            found_glmatl_ = find(glmatl_);            
            this.sv_ = {};
            [Nx,Ny,Nz] = this.registry.atlas_dims;
            glmmsk_3d = reshape(this.glmatl_, [Nx, Ny, Nz]); % in perceptron space, flip_{1,2}(NIfTI)
            R = this.registry.sphere_radius;
            I = this.registry.grid_spacing;
            isv = 1;
            
            for zr = 1:I:Nz
                for yr = 1:I:Ny
                    for xr = 1:I:Nx

                        if logical(glmmsk_3d(xr, yr, zr)) % assign spheres with centers inside glmmsk_3d                            
                            sphere_mask = zeros(Nx, Ny, Nz);
                            
                            for z = 1: Nz                        
                                for y = 1:Ny
                                    for x = 1:Nx
                                        
                                        if (x - xr)^2 + (y - yr)^2 + (z - zr)^2 < R^2
                                            sphere_mask(x,y,z) = 1; % set elements within ellipsoid to 1
                                        end
                                        
                                    end
                                end
                            end

                            sphere_mask = sphere_mask .* glmmsk_3d;
                            sphere_mask_glm = sphere_mask(found_glmatl_); %#ok<FNDSB> 
                            % numel(sphere_mask_glm) ~ 65549
                            % numel(sphere_mask) ~ 48*64*48
                            % numel(found_glmatl_) ~ 65549
                            % ordering of [x y z], [xr yr zr] is irrelevant
                            sphere_vox = find(sphere_mask_glm); % sphere_vox := indices
                            this.sv_{isv} = sphere_vox;
                            isv = isv+1;
                        end
                        
                    end
                end
            end
            
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

