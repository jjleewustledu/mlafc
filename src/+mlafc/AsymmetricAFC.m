classdef AsymmetricAFC < mlafc.AFC
	%% AsymmetricAFC  

	%  $Revision$
 	%  was created 04-Apr-2021 12:25:45 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlafc/src/+mlafc.
 	%% It was developed on Matlab 9.9.0.1592791 (R2020b) Update 5 for MACI64.  Copyright 2021 John Joowon Lee.
 	    
	properties (Dependent)
        BigBrain300
 		GMctx % ImagingContext2
        GMToYeo7_xfm
        N_BOLD
        sv
        Yeo7
        Yeo17
    end
    
    properties
        sl_fc_last
        tag
    end
    
    methods (Static)
        function brainNet_xi_x(coord, rsn, varargin)
            %% renders surfaces of fiber bundles for seed xi and base manifolds for voxel x.
            %  Requires Matlab R2014b.
            %  @param required coord in {'xi' 'x'}
            %  @param required i in {'smn' 'dmn' 'vis' 'fpc' 'lan' 'van' 'dan' }
            %  @param optional tag is char, e.g., '_noDelta'.
            
            ip = inputParser;
            addRequired(ip, 'coord', @ischar)
            addRequired(ip, 'rsn', @ischar)
            addOptional(ip, 'tag', '', @ischar)
            parse(ip, coord, rsn, varargin{:})
            ipr = ip.Results;
            
            workpath = pwd; % fullfile(getenv('HOME'), 'Box', 'Leuthardt_Epilepsy_Project', 'Segmentations_06_and_07_2020', '');
            
            if ~verLessThan('matlab', '8.5')
                error('mlafc:VersionError', 'AsymmetricAFC.brainNet_xi_x requires Matlab R2014')
            end
            assert(ischar(coord))
            assert(ischar(rsn))
            
            if strcmp(coord, 'xi')
                prefix = 'fiberbundle';
            else
                prefix = 'basemanifold';
            end
            BrainNet_MapCfg('~/MATLAB-Drive/BrainNetViewer_20191031/Data/SurfTemplate/BrainMesh_Ch2.nv', ...
                fullfile(workpath, sprintf('%s_%s_on_MNI152.node', prefix, rsn)), ...
                sprintf('%s_%s%s_on_MNI152.nii.gz', prefix, rsn, ipr.tag))
        end
        function brainNet_nifti(varargin)
            %% renders surfaces of fiber bundles for seed xi and base manifolds for voxel x.
            %  Requires Matlab R2014b.
            %  @param required fileprefix to NIfTI file on MNI152_1mm.
            %  @param required i in {'smn' 'dmn' 'vis' 'fpc' 'lan' 'van' 'dan' }
            %  @param optional tag is char, e.g., '_noDelta'.
            
            ip = inputParser;
            addRequired(ip, 'fileprefix', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            if ~verLessThan('matlab', '8.5')
                error('mlafc:VersionError', 'AsymmetricAFC.brainNet_xi_x requires Matlab R2014')
            end
            
            BrainNet_MapCfg('~/MATLAB-Drive/BrainNetViewer_20191031/Data/SurfTemplate/BrainMesh_Ch2.nv', ...
                sprintf('%s.nii.gz', ipr.fileprefix))
        end
        function buildNodeFiles(varargin)
            ip = inputParser;
            addParameter(ip, 'hemisphere', 'r', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            ipr.hemisphere = lower(ipr.hemisphere);
            
            rsns = {'smn' 'dmn' 'vis' 'fpc' 'van' 'dan' };
%             for r = 1:length(rsns)
%                 ic = mlfourd.ImagingContext2(sprintf('fiberbundle_%s_on_MNI152.nii.gz', rsns{r}));
%                 ic.fsleyes()
%             end
            
            % coords = [83 94 149;83 66 104; 83 47 77;44 150 99;54 150 71;35 94 84]; % R
            % coords = [99 91 151; 108 74 106; 115 51 74; 117 134 104; 54 150 71; 141 94 85]; % L
            if contains(ipr.hemisphere, 'r')
                coords1 = [-7 -32 77;-4 -60 32;-4 -79 5;-58 24 27;-55 22 -1;-29 -63 64];
            else
                coords1 = [9 -35 79; 18 -52 34;25 -75 2;27 8 32;-55 22 -1;51 -32 13];
            end
            x = -coords1(:,1);
            y = coords1(:,2);
            z = coords1(:,3);
            colors = ones(length(rsns), 1);
            size_xi = 12*ones(size(colors));
            size_x = 3*ones(size(colors));            
            labels = {'-' '-' '-' '-' '-' '-'};
            T_xi = table(x, y, z, colors, size_xi, labels');
            T_x = table(x, y, z, colors, size_x, labels');
            for r = 1:length(rsns)
                writetable(T_xi(r,:), ['fiberbundle_' rsns{r} '_on_MNI152.node'], 'FileType', 'text', 'WriteVariableNames', false, 'Delimiter', '\t')
                writetable(T_x(r,:), ['basemanifold_' rsns{r} '_on_MNI152.node'], 'FileType', 'text', 'WriteVariableNames', false, 'Delimiter', '\t')
            end
            
        end
        function buildFiberBundles()
            reg = mlafc.AFCRegistry.instance();
            sl_fc_mean_ic = mlafc.AsymmetricAFC.slfcMatToIC(reg.sl_fc_mean_mat, 'sl_fc_mean');
            fb = copy(sl_fc_mean_ic.nifti);

            fb_smn = copy(fb); fb_smn.img = fb.img(:,:,:,2433); fb_smn.fileprefix = 'fiberbundle_smn'; fb_smn.save
            fb_dmn = copy(fb); fb_dmn.img = fb.img(:,:,:,1835); fb_dmn.fileprefix = 'fiberbundle_dmn'; fb_dmn.save
            fb_vis = copy(fb); fb_vis.img = fb.img(:,:,:,1171); fb_vis.fileprefix = 'fiberbundle_vis'; fb_vis.save
            fb_fpc = copy(fb); fb_fpc.img = fb.img(:,:,:,1701); fb_fpc.fileprefix = 'fiberbundle_fpc'; fb_fpc.save
            fb_van = copy(fb); fb_van.img = fb.img(:,:,:,1007); fb_van.fileprefix = 'fiberbundle_van'; fb_van.save
            fb_dan = copy(fb); fb_dan.img = fb.img(:,:,:,2332); fb_dan.fileprefix = 'fiberbundle_dan'; fb_dan.save
            
            mlafc.AsymmetricAFC.flirt_to_MNI152('fiberbundle_smn.nii.gz')
            mlafc.AsymmetricAFC.flirt_to_MNI152('fiberbundle_dmn.nii.gz')
            mlafc.AsymmetricAFC.flirt_to_MNI152('fiberbundle_vis.nii.gz')
            mlafc.AsymmetricAFC.flirt_to_MNI152('fiberbundle_fpc.nii.gz')
            mlafc.AsymmetricAFC.flirt_to_MNI152('fiberbundle_van.nii.gz')
            mlafc.AsymmetricAFC.flirt_to_MNI152('fiberbundle_dan.nii.gz')
        end
        function buildFiberBundles_contrast(varargin)
            reg = mlafc.AFCRegistry.instance();
            
            ip = inputParser;
            addRequired(ip, 'matname', @isfile)
            addParameter(ip, 'sl_fc_mean_ic', [])
            addParameter(ip, 'hemisphere', 'r', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            if isempty(ipr.sl_fc_mean_ic)
                ipr.sl_fc_mean_ic = mlafc.AsymmetricAFC.slfcMatToIC(reg.sl_fc_mean_mat, 'sl_fc_mean');
            end
            
            ic = mlafc.AsymmetricAFC.slfcMatToIC(ipr.matname, 'sl_fc');
            ic = tanh(mlafc.AsymmetricAFC.atanh(ic) - mlafc.AsymmetricAFC.atanh(ipr.sl_fc_mean_ic));
            ic.filepath = pwd;
            fb = copy(ic.nifti);
            
            if contains(lower(ipr.hemisphere), 'r')
                fb_smn = copy(fb); fb_smn.img = fb.img(:,:,:,2433); fb_smn.fileprefix = 'fiberbundle_smn'; fb_smn.save
                fb_dmn = copy(fb); fb_dmn.img = fb.img(:,:,:,1835); fb_dmn.fileprefix = 'fiberbundle_dmn'; fb_dmn.save
                fb_vis = copy(fb); fb_vis.img = fb.img(:,:,:,1171); fb_vis.fileprefix = 'fiberbundle_vis'; fb_vis.save
                fb_fpc = copy(fb); fb_fpc.img = fb.img(:,:,:,1701); fb_fpc.fileprefix = 'fiberbundle_fpc'; fb_fpc.save
                fb_van = copy(fb); fb_van.img = fb.img(:,:,:,1007); fb_van.fileprefix = 'fiberbundle_van'; fb_van.save
                fb_dan = copy(fb); fb_dan.img = fb.img(:,:,:,2332); fb_dan.fileprefix = 'fiberbundle_dan'; fb_dan.save
            else
                fb_smn = copy(fb); fb_smn.img = fb.img(:,:,:,2392); fb_smn.fileprefix = 'fiberbundle_smn'; fb_smn.save
                fb_dmn = copy(fb); fb_dmn.img = fb.img(:,:,:,1618); fb_dmn.fileprefix = 'fiberbundle_dmn'; fb_dmn.save
                fb_vis = copy(fb); fb_vis.img = fb.img(:,:,:,1173); fb_vis.fileprefix = 'fiberbundle_vis'; fb_vis.save
                fb_fpc = copy(fb); fb_fpc.img = fb.img(:,:,:,1710); fb_fpc.fileprefix = 'fiberbundle_fpc'; fb_fpc.save
                fb_van = copy(fb); fb_van.img = fb.img(:,:,:,1007); fb_van.fileprefix = 'fiberbundle_van'; fb_van.save
                fb_dan = copy(fb); fb_dan.img = fb.img(:,:,:,2337); fb_dan.fileprefix = 'fiberbundle_dan'; fb_dan.save
            end
            
            mlafc.AsymmetricAFC.flirt_to_MNI152('fiberbundle_smn.nii.gz')
            mlafc.AsymmetricAFC.flirt_to_MNI152('fiberbundle_dmn.nii.gz')
            mlafc.AsymmetricAFC.flirt_to_MNI152('fiberbundle_vis.nii.gz')
            mlafc.AsymmetricAFC.flirt_to_MNI152('fiberbundle_fpc.nii.gz')
            mlafc.AsymmetricAFC.flirt_to_MNI152('fiberbundle_van.nii.gz')
            mlafc.AsymmetricAFC.flirt_to_MNI152('fiberbundle_dan.nii.gz')
        end
        function buildFiberBundles_nocontrast(varargin)
            ip = inputParser;
            addRequired(ip, 'matname', @isfile)
            addParameter(ip, 'hemisphere', 'r', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            ic = mlafc.AsymmetricAFC.slfcMatToIC(ipr.matname, 'sl_fc');
            ic.filepath = pwd;
            fb = copy(ic.nifti);
            
            if contains(lower(ipr.hemisphere), 'r')
                fb_smn = copy(fb); fb_smn.img = fb.img(:,:,:,2433); fb_smn.fileprefix = 'fiberbundle_smn_noDelta'; fb_smn.save
                fb_dmn = copy(fb); fb_dmn.img = fb.img(:,:,:,1835); fb_dmn.fileprefix = 'fiberbundle_dmn_noDelta'; fb_dmn.save
                fb_vis = copy(fb); fb_vis.img = fb.img(:,:,:,1171); fb_vis.fileprefix = 'fiberbundle_vis_noDelta'; fb_vis.save
                fb_fpc = copy(fb); fb_fpc.img = fb.img(:,:,:,1701); fb_fpc.fileprefix = 'fiberbundle_fpc_noDelta'; fb_fpc.save
                fb_van = copy(fb); fb_van.img = fb.img(:,:,:,1007); fb_van.fileprefix = 'fiberbundle_van_noDelta'; fb_van.save
                fb_dan = copy(fb); fb_dan.img = fb.img(:,:,:,2332); fb_dan.fileprefix = 'fiberbundle_dan_noDelta'; fb_dan.save
            else
                fb_smn = copy(fb); fb_smn.img = fb.img(:,:,:,2392); fb_smn.fileprefix = 'fiberbundle_smn_noDelta'; fb_smn.save
                fb_dmn = copy(fb); fb_dmn.img = fb.img(:,:,:,1618); fb_dmn.fileprefix = 'fiberbundle_dmn_noDelta'; fb_dmn.save
                fb_vis = copy(fb); fb_vis.img = fb.img(:,:,:,1173); fb_vis.fileprefix = 'fiberbundle_vis_noDelta'; fb_vis.save
                fb_fpc = copy(fb); fb_fpc.img = fb.img(:,:,:,1710); fb_fpc.fileprefix = 'fiberbundle_fpc_noDelta'; fb_fpc.save
                fb_van = copy(fb); fb_van.img = fb.img(:,:,:,1007); fb_van.fileprefix = 'fiberbundle_van_noDelta'; fb_van.save
                fb_dan = copy(fb); fb_dan.img = fb.img(:,:,:,2337); fb_dan.fileprefix = 'fiberbundle_dan_noDelta'; fb_dan.save
            end
            
            mlafc.AsymmetricAFC.flirt_to_MNI152('fiberbundle_smn_noDelta.nii.gz')
            mlafc.AsymmetricAFC.flirt_to_MNI152('fiberbundle_dmn_noDelta.nii.gz')
            mlafc.AsymmetricAFC.flirt_to_MNI152('fiberbundle_vis_noDelta.nii.gz')
            mlafc.AsymmetricAFC.flirt_to_MNI152('fiberbundle_fpc_noDelta.nii.gz')
            mlafc.AsymmetricAFC.flirt_to_MNI152('fiberbundle_van_noDelta.nii.gz')
            mlafc.AsymmetricAFC.flirt_to_MNI152('fiberbundle_dan_noDelta.nii.gz')
        end
        function buildPiFiberBundles(varargin)
            reg = mlafc.AFCRegistry.instance();
            
            ip = inputParser;
            addRequired(ip, 'matname', @isfile)
            addOptional(ip, 'sl_fc_mean', [], @isnumeric)
            addParameter(ip, 'pid', '', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            if isempty(ipr.sl_fc_mean)
                ld = load(reg.sl_fc_mean_mat, 'sl_fc_mean');
                ipr.sl_fc_mean = ld.sl_fc_mean;
            end
            
            ld = load(ipr.matname, 'sl_fc');
            sl_fc = ld.sl_fc;
            pi_energy = tanh(mean(mlafc.AsymmetricAFC.atanh(double(sl_fc)) - mlafc.AsymmetricAFC.atanh(double(ipr.sl_fc_mean)), 1, 'omitnan')); % 1 x 65546
            pi_energy_mat = sprintf('%s_pi_energy.mat', ipr.pid);
            save(pi_energy_mat, 'pi_energy')
            
            ic = mlafc.AsymmetricAFC.slfcMatToIC(pi_energy_mat, 'pi_energy');
            ifc = copy(ic.nifti); 
            ifc.fileprefix = 'pi_fiberbundles'; ifc.save            
            mlafc.AsymmetricAFC.flirt_to_MNI152('pi_fiberbundles.nii.gz')
        end
        function buildSoftmaxOnResection(varargin)
            ip = inputParser;
            addOptional(ip, 'toglob', 'PT*', @ischar)
            addParameter(ip, 'betFrac', 0.5, @isscalar)
            addParameter(ip, 'betZ', false, @islogical)
            addParameter(ip, 'cost', 'normmi', @ischar)
            addParameter(ip, 'dof', 12, @isscalar)
            addParameter(ip, 'search', 45, @isscalar)
            addParameter(ip, 'finesearch', 10, @isscalar)
            addParameter(ip, 'afcPrefix', 'AFC_product_tanhmae', @ischar)
            addParameter(ip, 'outTag', '_softmax', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            ipr.dof = num2str(ipr.dof);
            
            for g = globFoldersT(ipr.toglob)
                
                pwd0 = pushd(g{1});
                
                % ensure nii.gz
                if ~isfile([ipr.afcPrefix '.nii.gz'])
                    assert(isfile([ipr.afcPrefix '.4dfp.hdr']))
                    ic = mlfourd.ImagingContext2([ipr.afcPrefix '.4dfp.hdr']);
                    ic.nifti.save()
                end
                
                s = sprintf('-%i %i', ipr.search, ipr.search);
                fs = sprintf('%i', ipr.finesearch);
                opts = ['-bins 512 -cost ' ipr.cost ' -searchrx ' s ' -searchry ' s ' -searchrz ' s ' -dof ' ipr.dof ' -finesearch ' fs ' -interp trilinear'];                
                resection = globT([g{1} '_*_FLIRT_111.nii.gz']);
                assert(~isempty(resection))
                segmentation = globT([g{1} '_*_segmentation_final_111.nii.gz']);
                assert(~isempty(segmentation))
                
                r = '';
                try
                    % bet FLIRT_111
                    if ipr.betZ
                        resectionb = mlafc.AsymmetricAFC.betZ(resection{1}, ipr.betFrac);
                    else
                        resectionb = mlafc.AsymmetricAFC.bet(resection{1}, ipr.betFrac);
                    end
                    resectionb.fsleyes()
                    
                    % flirt AFC_product
                    [~,r] = mlbash(sprintf('flirt -in %s -ref %s -out %s%s -omat %s%s %s', ...
                        [ipr.afcPrefix '.nii.gz'], ...
                        resectionb.filename, ...
                        g{1}, [ipr.outTag '_111.nii.gz'], ...
                        g{1}, [ipr.outTag '_111.mat'], ...
                        opts));
                    
                    % quality control, viz.
                    mlafc.AsymmetricAFC.visualizeAfcProb(g{1})
                catch ME
                    warning('mlafc:RuntimeError', r)
                    handerror(ME)
                end
                
                popd(pwd0)
            end
        end    
        function Z_G = build_Z_G()
            %% cf. Bishop sec. 11.6.
            %  @return Z_G is numeric with size ~ 1 x Nvoxels
            
            import mlafc.AsymmetricAFC
            reg = mlafc.AFCRegistry.instance();
            sl_fc_gsp_ref_mat = @reg.sl_fc_gsp_ref_mat;
            betaT_ = reg.betaT;
            ldm = load(reg.sl_fc_mean_mat, 'sl_fc_mean');
            sl_fc_mean = ldm.sl_fc_mean;
            
            Z_G = zeros(size(ldm.sl_fc_mean));
            for r = 1:reg.ref_count
                %matfile = feval(sl_fc_gsp_ref_mat, r); %#ok<FVAL>
                ld = load(sl_fc_gsp_ref_mat(r), 'sl_fc_gsp_ref');
                energy = AsymmetricAFC.abs_Delta(ld.sl_fc_gsp_ref, sl_fc_mean);
                Z_G = Z_G + exp(-betaT_*energy);
            end
            
            save(reg.Z_G_mat, 'Z_G');
        end
        function garr = fullArrToMaskArr(farr)
            %% FULLARRTOMASKARR
            %  @param required farr is Nt x (48*64*48).
            %  @return garr is Nt x Nmask.
            
            assert(isnumeric(farr))
            assert(size(farr,2) == 48*64*48)
            assert(ismatrix(farr))
            
            mask_ = mlafc.AsymmetricAFC.readMaskArr();
            mask_(find(mask_)) = 1;   %#ok<FNDSB>
            found_mask_ = find(mask_); % Nmasked x 1
            
            Nt = size(farr,1);
            garr = zeros(Nt, length(found_mask_)); % Nt x Nmasked
            if 1 == Nt
                garr = farr(found_mask_');
                return
            end            
            for t = 1:Nt
                garr(t,:) = farr(t, found_mask_');
            end
        end 
        function farr = maskArrToFullArr(garr)
            %% MASKARRTOFULLARR
            %  @param required garr is Nt x 65549.
            %  @return farr is Nt x (48*64*48).
            
            assert(isnumeric(garr))
            assert(size(garr,2) == 65549 || size(garr,2) == 18611)
            assert(ismatrix(garr))        
            
            mask_ = mlafc.AsymmetricAFC.readMaskArr();
            mask_(find(mask_)) = 1;   %#ok<FNDSB>
            found_mask_ = find(mask_); % Nmasked x 1
            
            Nt = size(garr,1);
            farr = zeros(Nt, mlafc.AFCRegistry.instance().atlas_numel()); % Nt x 147456
            if 1 == Nt
                farr(found_mask_') = garr;
                return
            end            
            for t = 1:Nt
                farr(t,found_mask_') = garr(t,:);
            end
        end 
        function ic = maskArrToIC(varargin)
            ip = inputParser;
            addRequired(ip, 'arr', @(x) isvector(x))
            addParameter(ip, 'fileprefix', 'maskArrToIC', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;            
            
            arr = mlafc.AsymmetricAFC.maskArrToFullArr(ipr.arr);
            arr = arr';            
            img = reshape(arr, [48 64 48]);
            img(isnan(img)) = 0;
                         
            img = flip(img, 1);
            img = flip(img, 2);
            ifc = mlfourd.ImagingFormatContext(fullfile(getenv('REFDIR'), '711-2B_333.nii.gz'));
            ifc.img = img;
            ifc.filepath = pwd;
            ifc.fileprefix = ipr.fileprefix;
            ic = mlfourd.ImagingContext2(ifc);
        end
        function pi_fc = pi(fc)
            fcz = mlafc.AsymmetricAFC.atanh(double(fc)); % Fisher's z
            pi_fc = tanh(mean(fcz, 1, 'omitnan'));
        end
        function product = prepareMSC(varargin)
            %% Usage:  product createMscProbMaps(<subIdx>[, 'mscSubjectsMat', <mscSubjectsMat>])
            %  e.g.:   createMscProbMaps(6)
            %          createMscProbMaps(6, 'mscSubjectsMat', '/Users/jjlee/Box Sync/DeepNetFCProject/MSC/MSC_subject.mat')
            %
            %  @param subIdx is the subject index; or array of.
            %  @param mscSubjectsMat is a f.-q. filename; it must contain cell-array object MSC_subjects
            %         constructed as specified by mlnehorai.RyansMSC.createMscFourdfpMat().
            %  @return product is an mlfourd.ImagingFormatContext containing output from computeDenseConnectome();
            %          or array of.
            %  @return fullfile(mscSubjectPath, ['MSC_subject' num2str(subIdx) '_cat.mat']) containing
            %          object image := double(N_voxels, N_times), 
            %          with mscSubjectPath := fileparts(mscSubjectsMat); or array of.
            
            MSC_SUBJECTS_MAT = '/data/nil-bluearc/shimony/jjlee/MSC_Avis_denoising/MSC_subjects.mat';
            
            ip = inputParser;
            addRequired(ip, 'subIdx', @isnumeric)
            addParameter(ip, 'mscSubjectsMat', MSC_SUBJECTS_MAT, @isfile)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            mscDir = fileparts(ipr.mscSubjectsMat);
            pwd0 = pushd(mscDir);  
            load(ipr.mscSubjectsMat, 'MSC_subjects')
            
            product = cell(1, length(ipr.subIdx));
            for sidx = ipr.subIdx
                subPrefix = sprintf('MSC_subject%i', sidx);
                if ~isfile([subPrefix '_cat.mat'])
                    MSC_subject = MSC_subjects{sidx};
                    image = MSC_subject{1};
                    for ses = 2:length(MSC_subject)
                        image = cat(4, image, MSC_subject{ses});
                    end
                    image = reshape(image, [numel(image(:,:,:,1)) size(image,4)]);
                    save([subPrefix '_cat.mat'], 'image', '-v7.3')
                end            
%               product{sidx} = 
                mlperceptron.PerceptronFromMat.createProbMaps( ...
                    subPrefix, subPrefix, 'inmat', fullfile(mscDir, [subPrefix '_cat.mat']));  
            end
            if 1 == length(product)
                product = product{1};
            end
        
            popd(pwd0)
        end
        function arr = readMaskArr()
            reg = mlafc.AFCRegistry.instance();            
            if reg.GMonly
                arr = mlperceptron.PerceptronRegistry.read_N21_aparc_aseg_GMctx(); % 147456 x 1 =: 18611
            else
                arr = mlperceptron.PerceptronRegistry.read_glm_atlas_mask(); % 147456 x 1 =: 65549
            end            
        end
        function ic = slfcMatToIC(matname, objname, varargin)
            %% SLFCMATTOIC
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
            
            ld = load(matname, objname);
            obj = ld.(objname);
            assert(~isempty(obj))
            if 3 == ndims(obj)
                obj = squeeze(obj(ipr.refnum, :, :));
            end
            obj = mlafc.AsymmetricAFC.maskArrToFullArr(obj);
            obj = obj';
            
            Nspheres = size(obj, 2);
            if Nspheres > 1
                img = reshape(obj, [48 64 48 Nspheres]);
            else
                img = reshape(obj, [48 64 48]);
            end
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
        function visualizeAfcProb(varargin)
            ip = inputParser;
            addOptional(ip, 'toglob', 'PT*', @ischar)
            addParameter(ip, 'outTag', '_softmax', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            if strcmp(ipr.toglob, basename(pwd))
                visualize1(ipr.toglob)
                return
            end
            for g = globFoldersT(ipr.toglob)
                pwd0 = pushd(g{1});
                visualize1(g{1})
                popd(pwd0)
            end
            
            function visualize1(pid)
                resection = globT([pid '_*_FLIRT_111.nii.gz']);
                assert(~isempty(resection))
                segmentation = globT([pid '_*_segmentation_final_111.nii.gz']);
                assert(~isempty(segmentation))
                edge_seg = mlfourd.ImagingFormatContext(segmentation{1});
                edge_seg.img = edge3(edge_seg.img, 'approxcanny', 0.2);
                edge_seg.fileprefix = [strrep(segmentation{1}, '.nii.gz', '') '_edge'];
                if ~isfile(edge_seg.fqfilename)
                    edge_seg.save()
                end

                qc = mlfourd.ImagingContext2([pid ipr.outTag '_111.nii.gz']);
                qc.fsleyes(resection{1}, edge_seg.filename)
            end
        end
        function visualizeAfcThresh(varargin)
            ip = inputParser;
            addOptional(ip, 'toglob', 'PT*', @ischar)
            addOptional(ip, 'Nsigma', 2, @isscalar)
            addParameter(ip, 'outTag', '_tanhmae', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            if strcmp(ipr.toglob, basename(pwd))
                visualize1(ipr.toglob)
                return
            end
            for g = globFoldersT(ipr.toglob)
                pwd0 = pushd(g{1});
                visualize1(g{1})
                popd(pwd0)
            end
            
            function visualize1(pid)
                resection = globT([pid '_*_FLIRT_111.nii.gz']);
                assert(~isempty(resection))
                segmentation = globT([pid '_*_segmentation_final_111.nii.gz']);
                assert(~isempty(segmentation))
                edge_seg = mlfourd.ImagingFormatContext(segmentation{1});
                edge_seg.img = edge3(edge_seg.img, 'approxcanny', 0.2);
                edge_seg.fileprefix = [strrep(segmentation{1}, '.nii.gz', '') '_edge'];
                if ~isfile(edge_seg.fqfilename)
                    edge_seg.save()
                end

                product = mlfourd.ImagingContext2(sprintf('AFC_product%s.nii.gz', ipr.outTag));
                msk = product.binarized();
                img = product.nifti.img(logical(msk)); % 1 x Nmsk
                mu = mean(img);
                sigma = std(img);
                product = product.thresh(mu + ipr.Nsigma*sigma);
                product.fsleyes(resection{1}, edge_seg.filename)
            end
        end
        function x1 = x333_to_xMaskArr(varargin)
            %  @param x333 is [x,y,z] from the R^3 space of the 711-2B_333 atlas.
            %  @return x1 is the generalized coord in the R^1 space of this.glmatl_ or this.GMmsk_.
            
            ip = inputParser;
            addRequired(ip, 'x333', @(x) isvector(x))
            parse(ip, varargin{:})
            ipr = ip.Results;
            x = ipr.x333(1);
            y = ipr.x333(2);
            z = ipr.x333(3);
            
            img = zeros(48, 64, 48);
            img(x,y,z) = 1;
            img = flip(flip(img, 1), 2);
            farr = reshape(img, [1 48*64*48]);            
            garr = mlafc.AsymmetricAFC.fullArrToMaskArr(farr);
            x1 = find(garr);
            assert(isscalar(x1))
        end
    end

	methods        
        
        %% GET
        
        function g = get.BigBrain300(~)
            g = mlfourd.ImagingContext2( ...
                fullfile(getenv('REFDIR'), 'BigBrain300', 'BigBrain300_711-2b_allROIs.nii.gz'));
        end
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
        function g = get.N_BOLD(this)
            g = this.registry.N_BOLD;
        end
        function g = get.sv(this)
            g = this.sv_;
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
        
        function bold = applyMaskToBoldFrames(this, bold)
            if this.registry.GMonly
                mask = this.GMmsk_for_glm_;                
                bold = bold(:, logical(mask));
            end
        end
        function e = energy(this, fc)
            e = this.abs_Delta(fc, this.sl_fc_mean_);
        end
        function e = energy_similarity(this, fc)
            %% ENERGY_SIMILARITY <= 0.  Greater energy_similarity describes greater similarity between data and normal model.  
            
            switch this.registry.similarityTag
                case '_kldiv'
                    e = -this.kldiv(fc, this.sl_fc_mean_);
                case '_diag'
                    e = diag(this.acorrcoef(fc, this.sl_fc_mean_))';
                case '_tanhmae'
                    e = -tanh(mean(abs(this.atanh(double(fc)) - this.atanh(double(this.sl_fc_mean_))), 1, 'omitnan')); % mean abs error
                case '_mae'
                    e = -mean(abs(this.atanh(double(fc)) - this.atanh(double(this.sl_fc_mean_))), 1, 'omitnan'); % z-score
                case '_rmse'
                    e = -sqrt(mean((this.atanh(double(fc)) - this.atanh(double(this.sl_fc_mean_))).^2, 1, 'omitnan')); % z-score
                case '_mean2'
                    e = mean(this.atanh(double(this.acorrcoef(fc, this.sl_fc_mean_))), 2, 'omitnan')'; % z-score
                otherwise
                    error('mlafc:ValueError', 'AsymmetricAFC.energy_similarity')
            end
        end      
        function Z_E = Z_energy(this, fc)
            %% cf. Bishop sec. 11.6.
            %  @param required fc is functional connectivity, Gramian with size ~ Nseeds x Nvoxels
            %  @return Z_E_ is numeric with size ~ 1 x Nvoxels
            
            import mlafc.AsymmetricAFC
            reg = mlafc.AFCRegistry.instance();  
            sl_fc_gsp_ref_mat = @reg.sl_fc_gsp_ref_mat;
            betaT_ = reg.betaT;
            if ~isfile(reg.Z_G_mat)
                Z_G = this.build_Z_G();
            else
                ld = load(reg.Z_G_mat, 'Z_G');
                Z_G = ld.Z_G;
            end
            
            summand = zeros(size(this.sl_fc_mean_));
            for r = 1:reg.ref_count
                ld = load(sl_fc_gsp_ref_mat(r), 'sl_fc_gsp_ref');
                energy = AsymmetricAFC.abs_Delta(fc, ld.sl_fc_gsp_ref);
                summand = summand + exp(-betaT_*energy);
            end
                     
            % collect Z_E
            L = reg.ref_count;
            Z_E = summand .* Z_G / L;
            save(fullfile(this.patientdir, [this.patientid '_Z_E' reg.tag '.mat']), 'Z_E')
        end
        function this = explore_fc(this, varargin)
            %% EXPLORE_FC makes connectivity maps of all GSP subjects
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
                    bold_frames = this.applyMaskToBoldFrames(ld.bold_frames);

                    sl_fc_gsp_ref = zeros(Nsv__, this.N_BOLD); 
                    acorrcoef = @this.acorrcoef;
                    parfor isv = 1:Nsv__
                        bf_sv = mean(bold_frames(:, sv__{isv}), 2, 'omitnan'); %#ok<PFBNS> % N_t x 1
                        sl_fc_gsp_ref(isv,:) = ...
                            acorrcoef(bf_sv, bold_frames); % N_{sphere_vox} x 65549
                    end
                    save(reg.sl_fc_gsp_ref_mat(refnum), 'sl_fc_gsp_ref')
                    if reg.tanh_sandwich   
                        sl_fc_accum = sl_fc_accum + this.atanh(sl_fc_gsp_ref);
                    else
                        sl_fc_accum = sl_fc_accum + sl_fc_gsp_ref;
                    end
                    accum = accum + 1;
                catch ME
                    handwarning(ME)
                end
            end
            
            %% make mean map using GSP subjects
            
            if reg.tanh_sandwich  
                sl_fc_mean = tanh(sl_fc_accum/accum);
            else
                sl_fc_mean = sl_fc_accum/accum;
            end
            
            %% finalize          
            
            save(reg.sl_fc_mean_mat, 'sl_fc_mean')
            this.sl_fc_mean_ = sl_fc_mean;
            clear('sl_fc_mean')
            
            popd(pwd0)
        end  
        function this = explore_fc2(this, varargin)
            %% EXPLORE_FC2 makes connectivity maps of all GSP subjects using intermediates
            
            reg = this.registry;
            sv__ = this.sv;
            Nsv__ = length(sv__);
            
            sl_fc_accum = zeros(Nsv__, this.N_BOLD); 
            accum = 0;
            for refnum = 1:reg.ref_count
                try
                    ld = load(reg.sl_fc_gsp_ref_mat(refnum), 'sl_fc_gsp_ref');
                    sl_fc_gsp_ref = ld.sl_fc_gsp_ref;
                    if reg.tanh_sandwich                       
                        sl_fc_accum = sl_fc_accum + this.atanh(sl_fc_gsp_ref);
                    else
                        sl_fc_accum = sl_fc_accum + ld.sl_fc_gsp_ref;
                    end
                    accum = accum + 1;
                catch ME
                    handwarning(ME)
                end
            end
            
            %% make mean map using GSP subjects
            
            if reg.tanh_sandwich
                sl_fc_mean = tanh(sl_fc_accum/accum);
            else
                sl_fc_mean = sl_fc_accum/accum;
            end
            
            %% finalize          
            
            save(reg.sl_fc_mean_mat, 'sl_fc_mean')
            this.sl_fc_mean_ = sl_fc_mean;
            clear('sl_fc_mean')
        end  
        function [this,ipr] = makeSoftmax(this, varargin)
            %% MAKESOFTMAX of dissimilarity requires completion of explore_fc() which stores 
            %  this.sl_fc_gsp_, this.sl_fc_mean_.
            %  @param patientdir is a folder.
            %  @param patientid is char.
            %  @param Nframes is numeric.  If ~empty(Nframes), inferences uses only specified Nframes.
            %  @returns instance of mlafc.AsymmetricAFC.
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
            bold_frames = this.applyMaskToBoldFrames(ld.bold_frames); 
                          % N_t x N_vxls ~ 2307 x 18611 single
                          % N_t is variable across studies
            if ~isempty(ipr.Nframes)
                bold_frames = bold_frames(1:ipr.Nframes, :);                
                fprintf('makeSoftmax.bold_frames.size -> %s\n', mat2str(size(bold_frames)))
            end           
            
            %% per patient, build connectivity map
            
            sv__ = this.sv; 
            Nsv__ = length(sv__); % ~ 18611
            
            sl_fc = zeros(Nsv__, this.N_BOLD, class(bold_frames)); 
            acorrcoef = @this.acorrcoef;
            if size(sl_fc,1) == size(sl_fc,2)
                sl_fc = acorrcoef(bold_frames, bold_frames); % this.N_BOLD x this.N_BOLD ~ 18611 x 18611
            else
                parfor isv = 1:Nsv__ % 1:N_{sphere_vox}
                    bf_sv = mean(bold_frames(:, sv__{isv}), 2, 'omitnan'); %#ok<PFBNS> % N_t x 1
                    sl_fc(isv, :) = ...
                        acorrcoef(bf_sv, bold_frames); % N_{sphere_vox} x 65549
                end
            end
            
            %% project fiber bundles to base manifold
            
            this.sl_fc_last = sl_fc;
            prob = exp(-this.registry.betaT*this.energy(sl_fc)); % this.N_BOLD x this.N_BOLD

            %% assemble softmax
            
            prob = prob ./ this.Z_energy(sl_fc); % this.N_BOLD x this.N_BOLD
            map = this.pi(prob); % 1 x this.N_BOLD
            map = this.maskArrToFullArr(map); % 1 x 147456
            this.afc_map = map; % ease QA
            
            %% save            
            p = this.product('map', map);
            p.nifti.save();
            
            %% save more
            map1 = this.pi(sl_fc); % 1 x this.N_BOLD
            map1 = this.maskArrToFullArr(map1); % 1 x 147456
            p = this.product('fileprefix', [this.patientid '_pi_sl_fc'], 'map', map1);
            p.nifti.save();  
            save(fullfile(this.patientdir, [this.patientid '_sl_fc.mat']), 'sl_fc')
        end
        function ic = mapOfSpheres(this)
            reg = this.registry;
            glmarr = zeros(1, this.N_BOLD);
            for isv = 1:length(this.sv)
                glmarr(this.sv{isv}') = isv;
            end
            fullarr = this.maskArrToFullArr(glmarr);
            
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
        function p = product(this, varargin)
            ip = inputParser;
            addParameter(ip, 'fileprefix', [this.patientid '_softmax'], @ischar)
            addParameter(ip, 'map', this.afc_map, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            perreg = mlperceptron.PerceptronRegistry.instance();
            [d1,d2,d3] = perreg.atlas_dims;
            img = reshape(ipr.map, [d1 d2 d3]);
            img(isnan(img)) = 0;
            warning('mlafc:RuntimeWarning', ...
                'AFC.product() is flipping axes 1 & 2 to correct flips from mlperceptron.Fourdfp.Read4dfp')
            img = flip(flip(img, 1), 2);
            
            ifc = mlfourd.ImagingFormatContext( ...
                fullfile(getenv('REFDIR'), '711-2B_333.4dfp.hdr'));
            ifc.img = img;
            ifc.filepath = pwd;
            ifc.fileprefix = [ipr.fileprefix this.registry.fileTag];
            p = mlfourd.ImagingContext2(ifc);
        end
        function fc1 = resampleFunctionalConnectivity(this, varargin)
            %% resamples fc using samples{1,2} derived from mlpetersen.BigBrain300.imagingContext, typically
            %  the output of mlpetersen.BigBrain300.imagingContextSampling.  The contents of samples{1,2} determine
            %  the ordering by which fc1 has its indices sorted.  
            %  @param required fc is numeric for Nsv x Nbrain
            %  @param required samples1 is understood by ImagingContext2 and describes samples for xi.
            %  @param required samples2 is understood by ImagingContext2 and describes samples for x. 
            %  @returns fc1, which is smaller than fc, and which has indices sorted according to the
            %           contents of samples{1,2}.
            
            ip = inputParser;
            addRequired(ip, 'fc', @isnumeric)
            addRequired(ip, 'samples1', @(x) ~isempty(x))
            addRequired(ip, 'samples2', @(x) ~isempty(x))
            parse(ip, varargin{:})
            ipr = ip.Results;
            ipr.samples1 = mlfourd.ImagingContext2(ipr.samples1);
            ipr.samples2 = mlfourd.ImagingContext2(ipr.samples2);
            mask = mlafc.AsymmetricAFC.imagingContext_mask();
            
            % build indices1, xi
            samples1_img = ipr.samples1.fourdfp.img;
            samples1_img = flip(samples1_img, 1);
            samples1 = samples1_img(logical(mask)); % 1 x Nmask
            [sorted1_,indices1] = sort(samples1);
            indices1 = indices1(sorted1_ > 0);
            N1 = length(indices1);
            if size(ipr.fc, 1) == size(ipr.fc, 2)                
                fc_indices1_ = indices1';
            else
                fc_indices1_ = ones(1, N1);
                for i1 = 1:N1
                    for ixi = 1:length(this.sv)
                        roi_indices = unique(samples1(this.sv{ixi}));
                        if any(roi_indices == i1)
                            fc_indices1_(i1) = ixi;
                        end
                    end
                end
            end
            
            % build indices2, x         
            samples2_img = ipr.samples2.fourdfp.img;
            samples2_img = flip(samples2_img, 1);
            samples2 = samples2_img(logical(mask)); % 1 x Nmask  
            [sorted2_,indices2] = sort(samples2);
            indices2 = indices2(sorted2_ > 0);
            N2 = length(indices2);
            fc_indices2_ = indices2';
            
            % build fc1
            fc1 = zeros(N1, N2);
            for j = 1:N2
                for i = 1:N1
                    fc1(i, j) = ipr.fc(fc_indices1_(i), fc_indices2_(j));
                end
            end
        end
        function savefig(this, varargin)
            ip = inputParser;
            addRequired(ip, 'handle', @ishandle) % fig handle
            parse(ip, varargin{:})
            ipr = ip.Results;
                    
            dbs = dbstack;
            client_ = strrep(dbs(2).name, '.', '_');
            try
                savefig(ipr.handle, ...
                    fullfile(this.patientdir, sprintf('%s_%s%s.fig', this.patientid, client_, this.registry.tag)));
                figs = get(0, 'children');
                saveas(figs(1), ...
                    fullfile(this.patientdir, sprintf('%s_%s%s.png', this.patientid, client_, this.registry.tag)));
            catch ME
                handwarning(ME)
            end
        end
        function savefig_highres(this, varargin)
            ip = inputParser;
            addRequired(ip, 'handle', @ishandle) % fig handle
            addParameter(ip, 'position', [0 0 3500 2500], @isnumeric)
            addParameter(ip, 'res', 600, @isscalar) % fig res in dpi
            parse(ip, varargin{:})
            ipr = ip.Results;
            fh = ipr.handle;
            
            % set all units inside figure to normalized so that everything is scaling accordingly
            set(findall(fh,'Units','pixels'),'Units','normalized');
            
            % do not show figure on screen
            set(fh, 'visible', 'off')
            
            % set figure units to pixels & adjust figure size
            fh.Units = 'pixels';
            fh.OuterPosition = ipr.position;
            
            % recalculate figure size to be saved
            set(fh,'PaperPositionMode','manual')
            fh.PaperUnits = 'inches';
            fh.PaperPosition = ipr.position/ipr.res;
            
            % save figure
            dbs = dbstack;
            client_ = strrep(dbs(2).name, '.', '_');
            fileprefix = ...
                fullfile(this.patientdir, sprintf('%s_%s%s', this.patientid, client_, this.registry.tag));
            print(fh,fileprefix,'-dpng',sprintf('-r%d',ipr.res));
        end
        function fc1 = visualize_Deltaz_sl_fc(this)
            globbed = globT('*_sl_fc.mat');
            assert(~isempty(globbed))
            assert(~isempty(this.sl_fc_mean_))
            ld = load(globbed{1}, 'sl_fc');
            Deltaz = tanh(this.atanh(ld.sl_fc) - this.atanh(this.sl_fc_mean_));
            fc1 = this.visualize_sl_fc(Deltaz, 'tag', '_Deltaz_fc1');
        end
        function fc1 = visualize_sl_fc(this, varargin)
            ip = inputParser;
            addOptional(ip, 'sl_fc', this.sl_fc_last, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            bb300 = mlpetersen.BigBrain300('711-2B');
            s300 = bb300.imagingContextSampling(301);
            %s900 = bb300.imagingContextSampling(901);
            fc1 = this.resampleFunctionalConnectivity(ipr.sl_fc, s300, s300);
            
            % expensive to build, so save  
            save(sprintf('%s_visualize_sl_fc%s.mat', this.patientid, this.registry.tag), 'fc1')
            
            %surf(fc1)
            %hold on
            h = figure;
            imagesc(fc1)
            set(gca, 'FontSize', 32)
            colormap('jet')
            colorbar
            %hold off
            h1 = gca; 
            h1.XAxis.TickLength = [0 0];
            h1.YAxis.TickLength = [0 0];
            xlabel('voxels $\mathbf{x}$', 'Interpreter', 'latex', 'FontSize', 50)
            ylabel('seeds $\mathbf{\xi}$', 'Interpreter', 'latex', 'FontSize', 50)  
            %this.savefig_highres(h);
            this.savefig(h);
        end
        function pi_energy = visualize_pi_energy(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fc1', @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            h = figure;
            set(gca, 'FontSize', 32)
            pi_energy = tanh(mean(this.atanh(ipr.fc1), 1, 'omitnan'));
            max_ = max(pi_energy);
            min_ = min(pi_energy);
            range_ = max_ - min_;
            scaled_pi_energy = (pi_energy - min_)/range_;
            idx_pi_energy = uint8(255*scaled_pi_energy) + 1;
            jet256 = uint8(255*jet(256));
            cdat = jet256(idx_pi_energy, :);
            cdat1 = [cdat ones(300, 1, 'uint8')].';
            p = plot(1:300, pi_energy, 'r', 'LineWidth', 3);
            xlabel('voxels $\mathbf{x}$', 'Interpreter', 'latex', 'FontSize', 50)
            ylabel('$\mathbf{E}_{z,\xi} \Delta_z$ corr$(\mathbf{x}, \mathbf{\xi})$', 'Interpreter', 'latex', 'FontSize', 50)

            drawnow
            set(p.Edge, 'ColorBinding', 'interpolated', 'ColorData', cdat1)
            this.savefig(h);
        end
        
 		function this = AsymmetricAFC(varargin)
 			%% AsymmetricAFC
            %  @GMonly restricts to Carl's GM mask.
            %  @betaT is inverse temperature, default := 1.

 			this = this@mlafc.AFC(varargin{:});            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'patientid', '', @ischar)
            addParameter(ip, 'patientdir', pwd, @isfolder)
            addParameter(ip, 'similarityTag', '', @ischar)
            addParameter(ip, 'GMonly', true, @islogical)
            addParameter(ip, 'betaT', 0.5, @isscalar)
            addParameter(ip, 'tag', '_betaT0p5', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.registry_.min_num_vox = 1;
            this.registry_.tanh_sandwich = true;
            this.registry_.similarityTag = ipr.similarityTag;
            this.registry_.tag = ipr.tag;
            this.registry_.GMonly = ipr.GMonly;
            this.registry.betaT = ipr.betaT;
            this.patientid_ = ipr.patientid;
            this.patientdir = ipr.patientdir;
            
            this = sv_initialization(this);
            if isempty(this.sl_fc_mean_)
                this = load(this);
            end
        end        
    end 
    
    %% PROTECTED
    
    properties (Access = protected)        
        sv_ % 1 x 2825  cells with variable size [(x > 1) 1]
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
            reg = this.registry;
            if isfile(reg.sl_fc_mean_mat)
                mat = load(reg.sl_fc_mean_mat);
                this.sl_fc_mean_ = mat.sl_fc_mean;
            end
        end
        function this = sv_initialization(this)
            %% SV_INITIALIZATION
            %  @return this.sv_.
            
            % 04/18/19 KP
            % 06/06/19 JJL   
            
            %% Store sphere voxels indices
            
            mask_ = mlafc.AsymmetricAFC.readMaskArr();
            mask_(find(mask_)) = 1; %#ok<FNDSB>
            found_mask_ = find(mask_);            
            this.sv_ = {};
            [Nx,Ny,Nz] = this.registry.atlas_dims;
            mask_3d = reshape(mask_, [Nx, Ny, Nz]); % in perceptron space, flip_{1,2}(NIfTI)
            R = this.registry.sphere_radius;
            I = this.registry.grid_spacing;
            isv = 1;
            
            for zr = 1:I:Nz
                for yr = 1:I:Ny
                    for xr = 1:I:Nx

                        if logical(mask_3d(xr, yr, zr)) % assign spheres with centers inside glmmsk_3d                            
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

                            sphere_mask = sphere_mask .* mask_3d;
                            sphere_mask_glm = sphere_mask(found_mask_); %#ok<FNDSB> 
                            % numel(sphere_mask_glm) ~ 65549 | 18611
                            % numel(sphere_mask) ~ 48*64*48
                            % numel(found_mask_) ~ 65549 | 18611
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
    
    %% HIDDEN
  
    methods (Static, Hidden)
        function brainNet_xi_x_legacy(coord, i)
            %% renders surfaces of fiber bundles for seed xi and base manifolds for voxel x.
            %  Requires Matlab R2014b.
            %  @param required coord in {'xi' 'x'}
            %  @param required i in {'smn' 'dmn' 'vis' 'fpc' 'lan' 'van' 'dan' }
            
            if ~verLessThan('matlab', '8.5')
                error('mlafc:VersionError', 'AsymmetricAFC.brainNet_xi_x requires Matlab R2014')
            end
            assert(ischar(coord))
            assert(isscalar(i))
            
            xis = [1229 1463 1475 1695 1707 1909 1921 1933 2138 712];
            xpos = {[21 48 24] [21 48 27] [21 44 27] [21 45 30] [21 42 30] [21 42 33] [21 39 33] [21 36 33] [21 30 36] [21 4 15]};

            if strcmpi(coord, 'xi')
                BrainNet_MapCfg('~/MATLAB-Drive/BrainNetViewer_20191031/Data/SurfTemplate/BrainMesh_Ch2.nv', ...
                    sprintf('node%i_on_MNI152_xi.node', xis(i)), ...
                    sprintf('fiberbundle%i_on_MNI152.nii.gz', xis(i)))
            else
                x = xpos{i}(1) + 1;
                y = xpos{i}(2) + 1;
                z = xpos{i}(3) + 1;
                BrainNet_MapCfg('~/MATLAB-Drive/BrainNetViewer_20191031/Data/SurfTemplate/BrainMesh_Ch2.nv', ...
                    sprintf('node%i_on_MNI152_x.node', xis(i)), ...
                    sprintf('selectedBaseManifold_x_%i_%i_%i_on_MNI152.nii.gz', x, y, z))
            end
        end
        function buildBaseManifolds(varargin)
            %% BUILDBASEMANIFOLDS builds representations of base manifolds for visualization.
            %  BrainNet node center ~ [90 126 72]
            
            ip = inputParser;
            addParameter(ip, 'hemisphere', 'r', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            ipr.hemisphere = lower(ipr.hemisphere);            
            
            rsns = {'smn' 'dmn' 'vis' 'fpc' 'van' 'dan' };            
            if contains(ipr.hemisphere, 'r')
                x = {[27 23 44] [25 15 31] [26 8 21] [38 41 30] [35 41 22] [32 15 39]}; % fsleyes coords 
                % bb300_index = [44 118 161 199 215 265];
            else
                x = {[16 24 43] [21 16 28] [21 8 22] [10 42 30] [35 41 22] [18 15 40]}; % fsleyes coords
            end

            for r = 1:length(rsns)
                ic = mlafc.AsymmetricAFC.buildBaseManifold(rsns{r}, x{r} + 1);
                %ic.fsleyes
                ic.save
                fpout = strsplit(ic.fileprefix, '_x_');
                fpout = [fpout{1} '_on_MNI152'];
                mlafc.AsymmetricAFC.flirt_to_MNI152(ic.filename, 'fpout', fpout); %, 'interp', 'nearestneighbour');
            end
        end
        function ic = buildBaseManifold(varargin)
            %% builds base manifold for rsn using voxel coordinate x.
            %  @param required x is coordinate vector for 711-2B_333.
            %  @param radius is integer.
            %  @param stride is integer.
            %  @returns ImagingContext2 on 333-atlas with functional connectivity of x in base manifold
            %           to all fiber bundles defined by sv_radius*_stride*.mat
            
            ip = inputParser;
            addRequired(ip, 'rsn', @ischar)
            addRequired(ip, 'x', @isvector)
            addParameter(ip, 'radius', 2, @isscalar)
            addParameter(ip, 'stride', 3, @isscalar)
            addParameter(ip, 'ref_count', 500, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            ld = load(sprintf('sv_radius%i_stride%i.mat', ipr.radius, ipr.stride), 'sv');
            sv = ld.sv;
            ld = load(sprintf('sl_fc_mean_radius%i_stride%i_N%i_AsymmetricAFC.mat', ipr.radius, ipr.stride, ipr.ref_count), 'sl_fc_mean');
            sl_fc_mean = ld.sl_fc_mean;

            x1 = mlafc.AsymmetricAFC.x333_to_xMaskArr(ipr.x);
            sl_fc_xi = sl_fc_mean(:, x1)'; % as row
            garr = zeros(1, 65549);
            for isv = 1:length(sv)
                sv_set = sv{isv};
                for isv_set = 1:length(sv_set)
                    garr(sv_set(isv_set)) = sl_fc_xi(isv);
                end
            end
            x = ipr.x(1);
            y = ipr.x(2);
            z = ipr.x(3);
            ic = mlafc.AsymmetricAFC.glmmskArrToIC( ...
                garr, 'fileprefix', sprintf('basemanifold_%s_x_%i_%i_%i', ipr.rsn, x, y, z));
        end
        function buildBaseManifolds_legacy()
            %% BUILDBASEMANIFOLDS builds representations of base manifolds for visualization.
            %  xi ~ 1229 (DMN) 1463 (DMN) 1475 (DMN) 1695 1707 (FPN) 1909 1921 1933 (CON) 2138 (SMN) 712 (VIS)
            %  x on MNI152 ~
            %  {[21 48 24] [21 48 27] [21 44 27] [21 45 30] [21 42 30] [21 42 33] [21 39 33] [21 36 33] [21 30 36] [21 4 15]}
            %  BrainNet node center ~ [90 126 72]
            
            xpos = {[103 167 80] [103 167 90] [103 155 90] [103 157 100] [103 148 100] [103 148 110] [103 140 110] [103 130 110] [103 110 120] [103 31 55] ...
                    [20 42 36]}; 
                    % xi ~ 2091 (FPC), x on MNI ~ []

            for p = 1:length(xpos)
                ic = mlafc.AsymmetricAFC.buildBaseManifold_x(xpos{p} + 1);
                %ic.fsleyes
                ic.save
                mlafc.AsymmetricAFC.flirt_to_MNI152(ic.filename); %, 'interp', 'nearestneighbour');
            end
        end
        function ic = buildBaseManifold_x_legacy(varargin)
            %% builds base manifold selected by voxel coordinate x.
            %  @param required x is coordinate vector for 711-2B_333.
            %  @param radius is integer.
            %  @param stride is integer.
            %  @returns ImagingContext2 on 333-atlas with functional connectivity of x in base manifold
            %           to all fiber bundles defined by sv_radius*_stride*.mat
            
            ip = inputParser;
            addRequired(ip, 'x', @isvector)
            addParameter(ip, 'radius', 2, @isscalar)
            addParameter(ip, 'stride', 3, @isscalar)
            addParameter(ip, 'ref_count', 500, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            ld = load(sprintf('sv_radius%i_stride%i.mat', ipr.radius, ipr.stride), 'sv');
            sv = ld.sv;
            ld = load(sprintf('sl_fc_mean_radius%i_stride%i_N%i_AsymmetricAFC.mat', ipr.radius, ipr.stride, ipr.ref_count), 'sl_fc_mean');
            sl_fc_mean = ld.sl_fc_mean;

            x1 = mlafc.AsymmetricAFC.x333_to_xMaskArr(ipr.x);
            sl_fc_xi = sl_fc_mean(:, x1)'; % as row
            garr = zeros(1, 65549);
            for isv = 1:length(sv)
                sv_set = sv{isv};
                for isv_set = 1:length(sv_set)
                    garr(sv_set(isv_set)) = sl_fc_xi(isv);
                end
            end
            x = ipr.x(1);
            y = ipr.x(2);
            z = ipr.x(3);
            ic = mlafc.AsymmetricAFC.glmmskArrToIC( ...
                garr, 'fileprefix', sprintf('baseManifold_x_%i_%i_%i', x, y, z));
        end
        function buildFiberBundles_legacy()
            sl_fc_mean = mlfourd.ImagingContext2('sl_fc_mean_radius2_stride3_N500_AsymmetricAFC.nii.gz');
            fb = copy(sl_fc_mean.nifti);

            fb1229 = copy(fb); fb1229.img = fb.img(:,:,:,1229); fb1229.fileprefix = 'fiberbundle1229'; fb1229.save
            fb1463 = copy(fb); fb1463.img = fb.img(:,:,:,1463); fb1463.fileprefix = 'fiberbundle1463'; fb1463.save
            fb1475 = copy(fb); fb1475.img = fb.img(:,:,:,1475); fb1475.fileprefix = 'fiberbundle1475'; fb1475.save
            fb1695 = copy(fb); fb1695.img = fb.img(:,:,:,1695); fb1695.fileprefix = 'fiberbundle1695'; fb1695.save
            fb1707 = copy(fb); fb1707.img = fb.img(:,:,:,1707); fb1707.fileprefix = 'fiberbundle1707'; fb1707.save
            fb1909 = copy(fb); fb1909.img = fb.img(:,:,:,1909); fb1909.fileprefix = 'fiberbundle1909'; fb1909.save
            fb1921 = copy(fb); fb1921.img = fb.img(:,:,:,1921); fb1921.fileprefix = 'fiberbundle1921'; fb1921.save
            fb1933 = copy(fb); fb1933.img = fb.img(:,:,:,1933); fb1933.fileprefix = 'fiberbundle1933'; fb1933.save
            fb2138 = copy(fb); fb2138.img = fb.img(:,:,:,2138); fb2138.fileprefix = 'fiberbundle2138'; fb2138.save
            fb1975 = copy(fb); fb1975.img = fb.img(:,:,:,1975); fb1975.fileprefix = 'fiberbundle1975'; fb1975.save
            fb1990 = copy(fb); fb1990.img = fb.img(:,:,:,1990); fb1990.fileprefix = 'fiberbundle1990'; fb1990.save
            fb1822 = copy(fb); fb1822.img = fb.img(:,:,:,1822); fb1822.fileprefix = 'fiberbundle1822'; fb1822.save
            fb712 = copy(fb); fb712.img = fb.img(:,:,:,712); fb712.fileprefix = 'fiberbundle712'; fb712.save
            
            % FPC
            fb2091 = copy(fb); fb2091.img = fb.img(:,:,:,2091); fb2091.fileprefix = 'fiberbundle2091'; fb2091.save
        end
        function ic = glmmskArrToIC(varargin)
            ip = inputParser;
            addRequired(ip, 'arr', @(x) isvector(x))
            addParameter(ip, 'fileprefix', 'glmmskArrToIC', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;            
            
            arr = mlafc.AsymmetricAFC.maskArrToFullArr(ipr.arr);
            arr = arr';            
            img = reshape(arr, [48 64 48]);
            img(isnan(img)) = 0;
                         
            img = flip(img, 1);
            img = flip(img, 2);
            ifc = mlfourd.ImagingFormatContext(fullfile(getenv('REFDIR'), '711-2B_333.nii.gz'));
            ifc.img = img;
            ifc.filepath = pwd;
            ifc.fileprefix = ipr.fileprefix;
            ic = mlfourd.ImagingContext2(ifc);
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
    end
    
    methods (Hidden)
        function [this,ipr] = makeSoftmax_legacy(this, varargin)
            %% MAKESOFTMAX of dissimilarity requires completion of explore_fc() which stores 
            %  this.sl_fc_gsp_, this.sl_fc_mean_.
            %  @param patientdir is a folder.
            %  @param patientid is char.
            %  @param Nframes is numeric.  If ~empty(Nframes), inferences uses only specified Nframes.
            %  @returns instance of mlafc.AsymmetricAFC.
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
            bold_frames = this.applyMaskToBoldFrames(ld.bold_frames); 
                          % N_t x N_vxls ~ 2307 x 65549 single | 2307 x 18611 single
                          % N_t is variable across studies
            if ~isempty(ipr.Nframes)
                bold_frames = bold_frames(1:ipr.Nframes, :);                
                fprintf('makeSoftmax.bold_frames.size -> %s\n', mat2str(size(bold_frames)))
            end           
            
            %% per patient, build connectivity map
            
            sv__ = this.sv; 
            Nsv__ = length(sv__); % ~2441
            
            sl_fc = zeros(Nsv__, this.N_BOLD, class(bold_frames)); 
            acorrcoef = @this.acorrcoef;
            parfor isv = 1:Nsv__ % 1:N_{sphere_vox}
                bf_sv = mean(bold_frames(:, sv__{isv}), 2, 'omitnan'); %#ok<PFBNS> % N_t x 1
                sl_fc(isv, :) = ...
                    acorrcoef(bf_sv, bold_frames); % N_{sphere_vox} x 65549
            end
            
            %% project downsampled fibers to base manifold, using energy of similarity of fc to mean field fc,
            %% obtaning prob of dissimilarity.
            
            this.sl_fc_last = sl_fc;
            prob = exp(-this.registyr.betaT*this.energy_similarity(sl_fc)); % 1 x this.N_BOLD

            %% assemble softmax
            
            sum_prob = prob + this.sum_prob_refs(); % init with patient
            map = prob ./ sum_prob; % Boltzmann distribution
            map = this.maskArrToFullArr(map); % 1 x 147456
            this.afc_map = map; % ease QA
            
            %% save            
            save(fullfile(this.patientdir, [this.patientid '_sl_fc.mat']), 'sl_fc')
            p = this.product('map', map);
            p.nifti.save();
            
            %% save more
            map1 = this.maskArrToFullArr(tanh(mean(this.atanh(sl_fc), 1))); % 1 x 147456
            this.afc_map = map1;            
            p = this.product('fileprefix', [this.patientid '_pi_sl_fc'], 'map', map1);
            p.nifti.save();
            map2 = this.maskArrToFullArr(tanh(mean(this.atanh(sl_fc) - this.atanh(this.sl_fc_mean_), 1))); % 1 x 147456
            p = this.product('fileprefix', [this.patientid '_pi_sl_fc-sl_fc_mean'], 'map', map2);
            p.nifti.save();
            map2a = this.maskArrToFullArr(tanh(mean(abs(this.atanh(sl_fc) - this.atanh(this.sl_fc_mean_)), 1))); % 1 x 147456
            p = this.product('fileprefix', [this.patientid '_pi_abs_sl_fc-sl_fc_mean'], 'map', map2a);
            p.nifti.save();
            map3 = this.maskArrToFullArr(this.energy_similarity(sl_fc)); % 1 x 147456
            p = this.product('fileprefix', [this.patientid '_energy'], 'map', map3);
            p.nifti.save();  
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
                    sum_prob = sum_prob + exp(-this.energy_similarity(ld.sl_fc_gsp_ref));
                catch ME
                    handwarning(ME)
                end
            end
            save(this.registry.sl_fc_gsp_sum_prob_mat, 'sum_prob')
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

