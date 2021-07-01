classdef SymmetricAFC < mlafc.AFC
	%% SYMMETRICAFC  

	%  $Revision$
 	%  was created 04-Apr-2021 12:25:45 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlafc/src/+mlafc.
 	%% It was developed on Matlab 9.9.0.1592791 (R2020b) Update 5 for MACI64.  Copyright 2021 John Joowon Lee.
 	    
	properties (Dependent)
        BigBrain300
 		GMctx % ImagingContext2
        GMToYeo7_xfm
        N_BOLD
        Nframes
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
                error('mlafc:VersionError', 'SymmetricAFC.brainNet_xi_x requires Matlab R2014')
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
                error('mlafc:VersionError', 'SymmetricAFC.brainNet_xi_x requires Matlab R2014')
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
            sl_fc_mean_ic = mlafc.SymmetricAFC.slfcMatToIC(reg.sl_fc_mean_mat, 'sl_fc_mean');
            fb = copy(sl_fc_mean_ic.nifti);

            fb_smn = copy(fb); fb_smn.img = fb.img(:,:,:,2433); fb_smn.fileprefix = 'fiberbundle_smn'; fb_smn.save
            fb_dmn = copy(fb); fb_dmn.img = fb.img(:,:,:,1835); fb_dmn.fileprefix = 'fiberbundle_dmn'; fb_dmn.save
            fb_vis = copy(fb); fb_vis.img = fb.img(:,:,:,1171); fb_vis.fileprefix = 'fiberbundle_vis'; fb_vis.save
            fb_fpc = copy(fb); fb_fpc.img = fb.img(:,:,:,1701); fb_fpc.fileprefix = 'fiberbundle_fpc'; fb_fpc.save
            fb_van = copy(fb); fb_van.img = fb.img(:,:,:,1007); fb_van.fileprefix = 'fiberbundle_van'; fb_van.save
            fb_dan = copy(fb); fb_dan.img = fb.img(:,:,:,2332); fb_dan.fileprefix = 'fiberbundle_dan'; fb_dan.save
            
            mlafc.SymmetricAFC.flirt_to_MNI152('fiberbundle_smn.nii.gz')
            mlafc.SymmetricAFC.flirt_to_MNI152('fiberbundle_dmn.nii.gz')
            mlafc.SymmetricAFC.flirt_to_MNI152('fiberbundle_vis.nii.gz')
            mlafc.SymmetricAFC.flirt_to_MNI152('fiberbundle_fpc.nii.gz')
            mlafc.SymmetricAFC.flirt_to_MNI152('fiberbundle_van.nii.gz')
            mlafc.SymmetricAFC.flirt_to_MNI152('fiberbundle_dan.nii.gz')
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
                ipr.sl_fc_mean_ic = mlafc.SymmetricAFC.slfcMatToIC(reg.sl_fc_mean_mat, 'sl_fc_mean');
            end
            
            ic = mlafc.SymmetricAFC.slfcMatToIC(ipr.matname, 'sl_fc');
            ic = tanh(mlafc.SymmetricAFC.atanh(ic) - mlafc.SymmetricAFC.atanh(ipr.sl_fc_mean_ic));
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
            
            mlafc.SymmetricAFC.flirt_to_MNI152('fiberbundle_smn.nii.gz')
            mlafc.SymmetricAFC.flirt_to_MNI152('fiberbundle_dmn.nii.gz')
            mlafc.SymmetricAFC.flirt_to_MNI152('fiberbundle_vis.nii.gz')
            mlafc.SymmetricAFC.flirt_to_MNI152('fiberbundle_fpc.nii.gz')
            mlafc.SymmetricAFC.flirt_to_MNI152('fiberbundle_van.nii.gz')
            mlafc.SymmetricAFC.flirt_to_MNI152('fiberbundle_dan.nii.gz')
        end
        function buildFiberBundles_nocontrast(varargin)
            ip = inputParser;
            addRequired(ip, 'matname', @isfile)
            addParameter(ip, 'hemisphere', 'r', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            ic = mlafc.SymmetricAFC.slfcMatToIC(ipr.matname, 'sl_fc');
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
            
            mlafc.SymmetricAFC.flirt_to_MNI152('fiberbundle_smn_noDelta.nii.gz')
            mlafc.SymmetricAFC.flirt_to_MNI152('fiberbundle_dmn_noDelta.nii.gz')
            mlafc.SymmetricAFC.flirt_to_MNI152('fiberbundle_vis_noDelta.nii.gz')
            mlafc.SymmetricAFC.flirt_to_MNI152('fiberbundle_fpc_noDelta.nii.gz')
            mlafc.SymmetricAFC.flirt_to_MNI152('fiberbundle_van_noDelta.nii.gz')
            mlafc.SymmetricAFC.flirt_to_MNI152('fiberbundle_dan_noDelta.nii.gz')
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
            pi_energy = tanh(mean(mlafc.SymmetricAFC.atanh(double(sl_fc)) - mlafc.SymmetricAFC.atanh(double(ipr.sl_fc_mean)), 1, 'omitnan')); % 1 x 65546
            pi_energy_mat = sprintf('%s_pi_energy.mat', ipr.pid);
            save(pi_energy_mat, 'pi_energy')
            
            ic = mlafc.SymmetricAFC.slfcMatToIC(pi_energy_mat, 'pi_energy');
            ifc = copy(ic.nifti); 
            ifc.fileprefix = 'pi_fiberbundles'; ifc.save            
            mlafc.SymmetricAFC.flirt_to_MNI152('pi_fiberbundles.nii.gz')
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
                        resectionb = mlafc.SymmetricAFC.betZ(resection{1}, ipr.betFrac);
                    else
                        resectionb = mlafc.SymmetricAFC.bet(resection{1}, ipr.betFrac);
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
                    mlafc.SymmetricAFC.visualizeAfcProb(g{1})
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
            
            import mlafc.SymmetricAFC
            reg = mlafc.AFCRegistry.instance();
            sl_fc_gsp_ref_mat = @reg.sl_fc_gsp_ref_mat;
            betaT_ = reg.betaT;
            ldm = load(reg.sl_fc_mean_mat, 'sl_fc_mean');
            sl_fc_mean = ldm.sl_fc_mean;
            
            Z_G = zeros(size(ldm.sl_fc_mean));
            for r = 1:reg.ref_count
                %matfile = feval(sl_fc_gsp_ref_mat, r); %#ok<FVAL>
                ld = load(sl_fc_gsp_ref_mat(r), 'sl_fc_gsp_ref');
                energy = SymmetricAFC.abs_Delta(ld.sl_fc_gsp_ref, sl_fc_mean);
                Z_G = Z_G + exp(-betaT_*energy);
            end
            
            save(reg.Z_G_mat, 'Z_G');
        end
        function v = calc_violinplot(varargin)
            %  @param transform is a function handle, e.g., @probToKLSample.
            %  e.g.:
            %  mlafc.AFC.calc_violinplot('PT*', 'ylabel', '$-\log(p/p_0)$', 'interpreter', 'latex', 'transform', @mlafc.AFC.probToKLSample)
            
            ip = inputParser;
            addOptional(ip, 'toglob', 'PT*', @ischar)
            addParameter(ip, 'filestring', '', @ischar)
            addParameter(ip, 'outTag', '_softmax', @ischar)
            addParameter(ip, 'ylim', [], @isnumeric)
            addParameter(ip, 'ylabel', '$-\log(P(x)/P_0)$')
            addParameter(ip, 'interpreter', 'latex', @ischar)
            addParameter(ip, 'transform', @mlafc.AFC.probToKLSample)
            parse(ip, varargin{:})
            ipr = ip.Results;
            globbed = globFoldersT(ipr.toglob);
            reg = mlafc.AFCRegistry.instance();
            
            %mask = mlfourd.ImagingContext2(fullfile(getenv('REFDIR'), 'gm3d.nii.gz'));
            %mask = mask.binarized();
            mask = mlafc.AFC.imagingContext_mask();
            mask = mask.binarized();
            Nvxl = dipsum(mask);
            Npts = length(globbed);
            data = zeros(Nvxl, Npts, 'single');
            
            if strcmp(computer, 'MACI64')
                home = '/Users/jjlee/Box/Leuthardt_Epilepsy_Project/Segmentations_06_and_07_2020';
            else
                home = '/data/nil-bluearc/shimony/jjlee/FocalEpilepsy';
            end
            for ig = 1:Npts
                try
                    cd(fullfile(home, globbed{ig}));
                    if isempty(ipr.filestring)
                        afc_prob = mlfourd.ImagingContext2(sprintf('%s_softmax%s%s.nii.gz', ...
                            globbed{ig}, reg.fileTag, reg.productTag)); 
                    else
                        afc_prob = mlfourd.ImagingContext2(sprintf('%s_%s.nii.gz', ...
                            globbed{ig}, ipr.filestring));
                    end
                    img = afc_prob.nifti.img;
                    fprintf('gm sum -> %g\n', dipsum(img))
                    data(:,ig) = reshape(img(logical(mask)), [Nvxl 1]); 
                    if ~isempty(ipr.transform)
                        data(:,ig) = ipr.transform(data(:, ig));
                    end
                catch ME
                    handwarning(ME)
                end
            end
            cd(home)

            h = figure;
            if contains(ipr.toglob, 'PT')          
                %labels = {'PT15 (Ia)' 'PT26 (IV)' 'PT28 (Ia)' 'PT29 (IV)' 'PT34 (IIa)' 'PT35 (IV)' 'PT36 (Ia)'};  
                labels = {'1 (Ia)' '5 (IV)' '2 (Ia)' '6 (IV)' '4 (IIa)' '7 (IV)' '3 (Ia)'};  
                ordering = [1 3 7 5 2 4 6]; % pt28, pt15, pt36, pt34, pt26, pt29, p35
                labels = labels(ordering);
                data(:,:) = data(:,ordering);
            else                
                labels = cellfun(@(x) strrep(x, '_Ses1', ''), globbed, 'UniformOutput', false);
            end
            
            v = violinplot(data, labels, ...
                'ShowData', true, 'ShowNotches', false, 'ViolinAlpha', 0.03, ...
                'BoxColor', [.1 .1 .1], 'EdgeColor', [0.8 0.8 0.8], 'Bandwidth', 0.005);
            jet = colormap('jet');
            dv = 255/(length(v) - 1);
            for iv = 1:length(v)
                v(iv).BoxWidth = 0.02;
                idxjet = floor((iv-1)*dv) + 1;
                v(iv).ViolinColor = jet(idxjet,:);
            end
            if ~isempty(ipr.ylim)
                ylim(ipr.ylim)
            end
            ax = gca;
            ax.YRuler.Exponent = 0;
            %ytickformat('%.3f')

            if contains(globbed{1}, 'PT')
                cohort = 'pts';
                set(gca, 'fontsize', 26)
                ylabel(ipr.ylabel, 'Interpreter', ipr.interpreter, 'FontSize', 42)
                xlabel('patient ID (Engel class)', 'FontSize', 38)
            elseif contains(globbed{1}, 'Sub')
                cohort = 'subs';
                set(gca, 'fontsize', 26)
                xtickangle(90)
                ylabel(ipr.ylabel, 'Interpreter', ipr.interpreter, 'FontSize', 42)
                xlabel('subject ID', 'FontSize', 38)
            else
                cohort = 'entities';
                set(gca, 'fontsize', 14)
                ylabel(ipr.ylabel, 'Interpreter', ipr.interpreter, 'FontSize', 18)
                xlabel('entity ID', 'FontSize', 18)
            end
            savefig(h, ...
                sprintf('calc_violinplot%s%s_%i%s.fig', reg.fileTag, reg.productTag, Npts, cohort))
            figs = get(0, 'children');
            saveas(figs(1), ...
                sprintf('calc_violinplot%s%s_%i%s.png', reg.fileTag, reg.productTag, Npts, cohort))
%                close(figs(1))
        end 
        function c = client(n)
            assert(isscalar(n))
            dbs = dbstack;
            c = strrep(dbs(n+1).name, '.', '_');
        end
        function garr = fullArrToMaskArr(farr)
            %% FULLARRTOMASKARR
            %  @param required farr is Nt x (48*64*48).
            %  @return garr is Nt x Nmask.
            
            assert(isnumeric(farr))
            assert(size(farr,2) == 48*64*48)
            assert(ismatrix(farr))
            
            mask_ = mlafc.SymmetricAFC.readMaskArr();
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
            
            mask_ = mlafc.SymmetricAFC.readMaskArr();
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
            
            arr = mlafc.SymmetricAFC.maskArrToFullArr(ipr.arr);
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
            fcz = mlafc.SymmetricAFC.atanh(double(fc)); % Fisher's z
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
            obj = mlafc.SymmetricAFC.maskArrToFullArr(obj);
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
        function visualizePT(varargin)
            ip = inputParser;
            addRequired(ip, 'pid', @(x) contains(x, 'PT'))
            addOptional(ip, 'Z_G', [], @isnumeric)
            addParameter(ip, 'KLSample', false, @islogical)
            parse(ip, varargin{:})
            ipr = ip.Results;
            if isempty(ipr.Z_G)
                ld = load('Z_G_radius1_stride1_N500_GMonly.mat', 'Z_G');
                ipr.Z_G = ld.Z_G;
            end
                        
            reg = mlafc.AFCRegistry.instance();
            pwd0 = pushd(fullfile(reg.workpath, ipr.pid));
            
            ld = load([ipr.pid '_sl_fc.mat'], 'sl_fc');
            jafc = mlafc.SymmetricAFC();
            if ipr.KLSample
                prob = exp(-jafc.energy(ld.sl_fc)) ./ ipr.Z_G;
                prob = jafc.probToKLSample(prob);
                jafc.visualize_sl_fc(prob, ...
                    'fileprefix', 'Gramian_KLSample', 'z', false, 'flip_colormap', false, 'caxis', [-.15 .2]); % [1.8e-3 2.3e-3]);
            else
                jafc.visualize_sl_fc(jafc.abs_Delta(ld.sl_fc, jafc.sl_fc_mean_), ...
                    'fileprefix', 'Gramian_abs_Delta_sl_fc_z', 'caxis', [0 0.5]);
            end
            
            popd(pwd0)
        end
        function visualizeKLSampleFsleyes(varargin)
            ip = inputParser;
            addOptional(ip, 'toglob', 'PT*', @ischar)
            addParameter(ip, 'outTag', '_softmax', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            reg = mlafc.AFCRegistry.instance();
            
            for g = globFoldersT(ipr.toglob)
                pwd0 = pushd(fullfile(reg.workpath, g{1}));
                in = sprintf('%s_softmax%s%s.nii.gz', g{1}, reg.fileTag, reg.productTag);
                assert(isfile(in))
                mat = sprintf('%s_softmax_111.mat', g{1});
                out = sprintf('%s_111', myfileprefix(in));
                ref = fullfile(getenv('REFDIR'), '711-2B_111.nii.gz');
                mlafc.AFC.applyxfm(in, mat, 'out_fqfp', out, 'icref', ref, 'interp', 'nearestneighbour') 
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

                bb300 = fullfile(getenv('REFDIR'), 'BigBrain300', 'BigBrain300_711-2b_allROIs_111.nii.gz');
                
                qc = mlfourd.ImagingFormatContext([out '.nii.gz']);
                img = qc.img(qc.img ~= 0);
                qc.img(qc.img ~= 0) = -log(500*img);
                qc.fileprefix = [qc.fileprefix '_logpp0'];
                qc.save()
                qc.fsleyes(resection{1}, edge_seg.filename, bb300)
            end
        end
        function visualizeRSNROIInSegmentation()
            import mlfourd.ImagingContext2
            reg = mlafc.AFCRegistry.instance();
            bb = ImagingContext2(fullfile(getenv('REFDIR'), 'BigBrain300', 'BigBrain300_711-2b_allROIs_111.nii.gz'));
            ld = load(fullfile(reg.workpath, 'inv_uhits.mat'));
            inv_uhits = ld.inv_uhits;

            globbed = globFoldersT('PT*');
            pids = [2 5 1 6 4 7 3];
            engels = {'Ia' 'IV' 'Ia' 'IV' 'IIa' 'IV' 'Ia'};
            for ig = 1:length(globbed)
                pwd0 = pushd(fullfile(reg.workpath, globbed{ig}));
                gm = ImagingContext2(sprintf('%s_softmax%s%s_111.nii.gz', globbed{ig}, reg.fileTag, reg.productTag));
                gm = gm.binarized();
                globbed_seg = glob([globbed{ig} '*_segmentation_final_111.nii.gz']);
                seg = ImagingContext2(globbed_seg{1});
                seg = seg.binarized();
                sample_ = seg .* gm .* bb;
                sample = sample_.nifti.img;
                sample1 = sample(sample > 0);
                sample2 = inv_uhits(sample1);
                figure
                histogram(sample2, 201)
                xlim([0 201])
                set(gca, 'FontSize', 32)
                set(gca, 'XTickLabel', [])
                %ylabel('Volume of RSN ROI within Resection (\muL)', 'FontSize', 50)
                %title(sprintf('patient %i (Engel %s)', pids(ig), engels{ig}), 'FontSize', 50)                
                h1 = gca;
                h1.XAxis.TickLength = [0 0];
                popd(pwd0)
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
            garr = mlafc.SymmetricAFC.fullArrToMaskArr(farr);
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
        function g = get.Nframes(this)
            g = this.Nframes_;
        end
        function g = get.N_BOLD(this)
            g = this.registry.N_BOLD;
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
                    error('mlafc:ValueError', 'SymmetricAFC.energy_similarity')
            end
        end      
        function Z_E = Z_energy(this, fc)
            %% cf. Bishop sec. 11.6.
            %  @param required fc is functional connectivity, Gramian with size ~ Nseeds x Nvoxels
            %  @return Z_E_ is numeric with size ~ 1 x Nvoxels
            
            import mlafc.SymmetricAFC
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
                energy = SymmetricAFC.abs_Delta(fc, ld.sl_fc_gsp_ref);
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
                     
            sl_fc_accum = zeros(this.N_BOLD, this.N_BOLD); 
            accum = 0;
            for refnum = 1:reg.ref_count
                try
                    refid = reg.gtm500_ids{refnum};
                    ld = load(fullfile(gtm500_dir, ...
                               refid, ...
                               'Perceptron', ...
                               sprintf('%s%s', refid, ref_resid_mat)), 'bold_frames'); % bold_frames has N_t x 18611  
                    bold_frames = this.applyMaskToBoldFrames(ld.bold_frames);
                    sl_fc_gsp_ref = this.acorrcoef(bold_frames, bold_frames); % 18611 x 18611
                        
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
            %% EXPLORE_FC2 makes a std map of all GSP subjects using existing intermediates
            
            reg = this.registry;
            
            sl_fc_accum = zeros(this.N_BOLD, this.N_BOLD); 
            accum = 0;
            for refnum = 1:reg.ref_count
                try
                    ld = load(reg.sl_fc_gsp_ref_mat(refnum), 'sl_fc_gsp_ref');
                    sl_fc_gsp_ref = ld.sl_fc_gsp_ref;
                    if reg.tanh_sandwich                       
                        sl_fc_accum = sl_fc_accum + ...
                            (this.atanh(sl_fc_gsp_ref) - this.atanh(this.sl_fc_mean_)).^2;
                    else
                        sl_fc_accum = sl_fc_accum + ...
                            (ld.sl_fc_gsp_ref - this.sl_fc_mean_).^2;
                    end
                    accum = accum + 1;
                catch ME
                    handwarning(ME)
                end
            end
            
            %% make mean map using GSP subjects
            
            if reg.tanh_sandwich
                sl_fc_std = tanh(sqrt(sl_fc_accum/(accum - 1)));
            else
                sl_fc_std = sqrt(sl_fc_accum/(accum - 1));
            end
            
            %% finalize          
            
            save(reg.sl_fc_std_mat, 'sl_fc_std')
            clear('sl_fc_std')
        end  
        function fp = fileprefix(this, varargin)
            % @param n is the index of the client in the dbstack; default := 1
            
            ip = inputParser;
            addParameter(ip, 'n', 1, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            fp = sprintf('%s_%s%s', this.patientid, this.client(ipr.n+1), this.registry.tag);
        end
        function fqfp = fqfileprefix(this, varargin)
            fqfp = fullfile(this.patientdir,this.fileprefix(varargin{:}));
        end
        function [this,ipr] = makeSoftmax(this, varargin)
            %% MAKESOFTMAX of dissimilarity requires completion of explore_fc() which stores 
            %  this.sl_fc_gsp_, this.sl_fc_mean_.
            %  @param patientdir is a folder.
            %  @param patientid is char.
            %  @param Nframes is numeric.  If ~empty(Nframes), inferences uses only specified Nframes.
            %  @returns instance of mlafc.SymmetricAFC.
            %  @returns inputParser.Results.
            %  @returns softmax in this.product as ImagingContext2.
            
            import mlpark.*
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'patientdir', @isfolder)
            addRequired( ip, 'patientid', @ischar)
            addParameter(ip, 'Nframes', this.Nframes, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;   
            this.patientdir = ipr.patientdir;
            this.patientid = ipr.patientid;
            
            %% per patients, use results of mlperceptron.PerceptronRelease factory
            
            mlperceptron.PerceptronFromFormat.createProbMaps( ...
                this.patientdir, this.patientid, 'frameCensoring', false)
            try
                ld = load( ...
                    fullfile(this.patientdir, 'Perceptron', ...
                    sprintf('%s%s', this.patientid, this.registry.perceptron_resid_mat)), ...
                    'bold_frames');
            catch ME
                handerror(ME)
                %ld = load( ...
                %    fullfile(this.patientdir, 'Perceptron', ...
                %    sprintf('%s%s', this.patientid, this.registry.perceptron_uout_resid_mat)), ...
                %    'bold_frames');
            end
            bold_frames = this.applyMaskToBoldFrames(ld.bold_frames); 
                          % N_t x N_vxls ~ 2307 x 18611 single
                          % N_t is variable across studies
            if ~isempty(ipr.Nframes)
                bold_frames = bold_frames(1:ipr.Nframes, :);                
                fprintf('makeSoftmax.bold_frames.size -> %s\n', mat2str(size(bold_frames)))
            end           
            
            %% per patient, build connectivity map
            
            acorrcoef = @this.acorrcoef;
            sl_fc = acorrcoef(bold_frames, bold_frames); % this.N_BOLD x this.N_BOLD ~ 18611 x 18611
            
            %% project fiber bundles to base manifold
            
            this.sl_fc_last = sl_fc;
            prob = exp(-this.registry.betaT*this.energy(sl_fc)); % this.N_BOLD x this.N_BOLD

            %% assemble softmax
            
            ld = load(this.registry.Z_G_mat, 'Z_G');
            prob = prob ./ ld.Z_G; % this.N_BOLD x this.N_BOLD
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
            sl_fc_file = fullfile(this.patientdir, [this.patientid '_sl_fc.mat']);
            %if ~isfile(sl_fc_file)
            save(sl_fc_file, 'sl_fc')
            %end
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
            ifc.fileprefix = [ipr.fileprefix this.registry.fileTag this.registry.productTag];
            p = mlfourd.ImagingContext2(ifc);
        end
        function fc1 = resampleFunctionalConnectivity(~, varargin)
            %% resamples fc using samples{1,2} derived from mlpetersen.BigBrain300.imagingContext, typically
            %  the output of mlpetersen.BigBrain300.imagingContextSampling.  The contents of samples{1,2} determine
            %  the ordering by which fc1 has its indices sorted.  
            %  @param required fc is numeric for Nvoxels x Nvoxels
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
            mask = mlafc.SymmetricAFC.imagingContext_mask();
            
            % build indices1, xi
            samples1_img = ipr.samples1.fourdfp.img;
            samples1_img = flip(samples1_img, 1);
            samples1 = samples1_img(logical(mask)); % 1 x Nmask
            [sorted1_,indices1] = sort(samples1);
            indices1 = indices1(sorted1_ > 0);
            N1 = length(indices1);               
            fc_indices1_ = indices1';
            
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
        function fc1 = resampleFunctionalConnectivityAll(this, varargin)
            %% resamples fc using samples{1,2} derived from mlpetersen.BigBrain300.imagingContext, typically
            %  the output of mlpetersen.BigBrain300.imagingContextSampling.  The contents of samples1 determine
            %  the ordering by which fc1 has its indices sorted.  
            %  @param required fc is numeric for Nvoxels x Nvoxels, in Carl Hacker's format.
            %  @param required samples1 is understood by ImagingContext2 and describes samples for x.
            %  @returns fc1, which is smaller than fc, and which has indices sorted according to the
            %           contents of samples1.
            
            ip = inputParser;
            addRequired(ip, 'fc', @isnumeric)
            addRequired(ip, 'samples1', @(x) ~isempty(x))
            parse(ip, varargin{:})
            ipr = ip.Results;
            ipr.samples1 = mlfourd.ImagingContext2(ipr.samples1);
            ipr.samples1 = flip(ipr.samples1, 1);
            ipr.samples1 = flip(ipr.samples1, 2); % nifti
            mask = mlafc.SymmetricAFC.imagingContext_mask();
            mask = flip(mask, 1); % 4dfp
            
            % build indices1, x
            samples1 = ipr.samples1.nifti.img(logical(mask)); % 1 x Nmask, containing [0:300]
            
            % build fc1
            ld = load(fullfile(this.registry.workpath, 'uhits.mat'), 'uhits');
            uhits = ld.uhits';
            Nhits = length(uhits);
            fc1 = zeros(Nhits, Nhits);
            for j = 1:Nhits
                for i = 1:Nhits
                    rows = find(samples1 == uhits(i));
                    cols = find(samples1 == uhits(j));
                    fc1(i,j) = tanh(mean(mean(atanh(ipr.fc(rows,cols)), 'omitnan'), 'omitnan'));
                end
            end
            
            fc1 = this.reshape201(fc1);
        end
        function fc1 = resampleFunctionalConnectivityReals(this, varargin)
            %% resamples fc using samples{1,2} derived from mlpetersen.BigBrain300.imagingContext, typically
            %  the output of mlpetersen.BigBrain300.imagingContextSampling.  The contents of samples1 determine
            %  the ordering by which fc1 has its indices sorted.  
            %  @param required fc is numeric for Nvoxels x Nvoxels, in Carl Hacker's format.
            %  @param required samples1 is understood by ImagingContext2 and describes samples for x.
            %  @returns fc1, which is smaller than fc, and which has indices sorted according to the
            %           contents of samples1.
            
            ip = inputParser;
            addRequired(ip, 'fc', @isnumeric)
            addRequired(ip, 'samples1', @(x) ~isempty(x))
            parse(ip, varargin{:})
            ipr = ip.Results;
            ipr.samples1 = mlfourd.ImagingContext2(ipr.samples1);
            ipr.samples1 = flip(ipr.samples1, 1);
            ipr.samples1 = flip(ipr.samples1, 2); % nifti
            mask = mlafc.SymmetricAFC.imagingContext_mask();
            mask = flip(mask, 1); % 4dfp
            
            % build indices1, x
            samples1 = ipr.samples1.nifti.img(logical(mask)); % 1 x Nmask, containing [0:300]
            
            % build fc1
            ld = load(fullfile(this.registry.workpath, 'uhits.mat'), 'uhits');
            uhits = ld.uhits';
            Nhits = length(uhits);
            fc1 = zeros(Nhits, Nhits);
            for j = 1:Nhits
                for i = 1:Nhits
                    rows = find(samples1 == uhits(i));
                    cols = find(samples1 == uhits(j));
                    fc1(i,j) = mean(mean(ipr.fc(rows,cols), 'omitnan'), 'omitnan');
                end
            end
            
            fc1 = this.reshape201(fc1);
        end
        function fc1 = reshape201(~, fc)
            %% arranges Gramian fc (201 x 201) to move DAN and inserts rows/cols nans to separate 15 RSNs
            
            assert(size(fc,1) == 201)
            assert(size(fc,2) == 201)
            
            % relocate dan
            ortho = eye(201);
            ortho1 = ortho;
            ortho1(180:196,180:196) = zeros(17); % cut
            ortho1(153:169,180:196) = ortho(180:196,180:196); % paste
            ortho1(153:179,153:179) = zeros(27);
            ortho1(170:196,153:179) = ortho(153:179,153:179);
            fc = ortho1 * fc * ortho1';
            
            % relocate rew
            ortho = eye(201);
            ortho1 = zeros(201);
            ortho1(197:201,1:5) = ortho(1:5,1:5);
            ortho1(1:196,6:201) = ortho(6:201,6:201);
            fc = ortho1 * fc * ortho1';
            
            % expand fc to show RSN boundary markers in black (nan)
            fc_ = zeros(215);
            fc_(1:201,1:201) = fc;            
            
            % define RSN ranges
            rr = containers.Map('KeyType', 'double', 'ValueType', 'any'); % RSN ranges
            rr(1) = 1:22;     % smd
            rr(2) = 23:26;    % smi
            rr(3) = 27:45;    % con
            rr(4) = 46:57;    % aud *
            rr(5) = 58:96;    % dmn *
            rr(6) = 97:99;    % pmn
            rr(7) = 100:120;  % vis *
            rr(8) = 121:134;  % fpn * 
            rr(9) = 135:142;  % sal
            rr(10) = 143:147; % van
            rr(11) = 148:164; % dan
            rr(12) = 165:179; % bga
            rr(13) = 180:191; % tha
            rr(14) = 192:196; % mtl
            rr(15) = 197:201; % rew
            
            % relocate RSNs and insert RSN boundary markers
            ortho = eye(215);
            ortho1 = zeros(215,215);
            for i = 0:14
                j = i + 1;
                rr_ = rr(j);
                ortho1(rr_+i,rr_) = ortho(rr_,rr_);
            end            
            % e.g., ortho1((1:22)+0,1:22) = ortho(1:22,1:22);
            %       ortho1((23:26)+1,23:26) = ortho(23:26,23:26);
            %       ortho1((27:45)+2,27:45) = ortho(27:45,27:45);
            fc1 = ortho1 * fc_ * ortho1'; 
            for j = 1:14
                rr_ = rr(j);
                fc1(rr_(end)+j,:) = nan(1,215);
                fc1(:,rr_(end)+j) = nan(215,1);
            end
        end
        function savefig(this, varargin)
            
            dbs = dbstack;
            client_ = strrep(dbs(2).name, '.', '_');
            fp = fullfile(this.patientdir, sprintf('%s_%s%s', this.patientid, client_, this.registry.tag));
            
            ip = inputParser;
            addRequired(ip, 'handle', @ishandle) % fig handle
            addParameter(ip, 'fileprefix', '', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            if isempty(ipr.fileprefix)
                ipr.fileprefix = fp;
            end
                    
            try
                savefig(ipr.handle, [ipr.fileprefix '.fig']);
                figs = get(0, 'children');
                saveas(figs(1), [ipr.fileprefix '.png']);
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
        function fc1 = visualize_Deltaz_sl_fc(this, varargin)
            globbed = globT('*_sl_fc.mat');
            assert(~isempty(globbed))
            assert(~isempty(this.sl_fc_mean_))
            ld = load(globbed{1}, 'sl_fc');
            Deltaz = tanh(this.atanh(ld.sl_fc) - this.atanh(this.sl_fc_mean_));            
            fc1 = this.visualize_sl_fc(Deltaz, 'n', 2, varargin{:});
        end
        function fc1 = visualize_sl_fc(this, varargin)
            %% visualizes Gramian for Carl's gray-matter mask.
            % @param n is the index of the client in the dbstack; default := 1; used with fileprefix(this).
            %
            % uhit coord	RSN (15)          ortho   ortho
            % 
            % 6:27		    smd                         1:22 
            % 28:31		    smi                        23:26
            % 32:50		    con                        27:45
            % 51:61		    aud                        46:56
            % 62:101		dmn                        57:96
            % 102:104		pmn                        97:99
            % 105:124		vis                       100:119
            % 125:139		fpn                       120:134
            % 140:147		sal                       135:142
            % 148:152		van                       143:147
            % 180:196		dan *             153:169 148:164
            % 153:167		bga               170:184 165:179
            % 168:179		tha               185:196 180:191
            % 197:198		amygdala    } MTL 197:201 192:196
            % 199:201		hippocampus }
            % 1:5			reward *                  197:201
            % * indicates relocation

            ip = inputParser;
            addOptional(ip, 'sl_fc', this.sl_fc_last, @isnumeric)
            addParameter(ip, 'fileprefix', '', @ischar)
            addParameter(ip, 'z', true, @islogical)
            addParameter(ip, 'caxis', [], @isnumeric)
            addParameter(ip, 'colormap', 'jet', @ischar)
            addParameter(ip, 'flip_colormap', false, @islogical)
            addParameter(ip, 'save', true, @islogical)
            addParameter(ip, 'n', 1, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            bb300 = mlpetersen.BigBrain300('711-2B');
            %s300 = bb300.imagingContextSampling(301);
            sall = bb300.imagingContextSampling(20*300);
            fc1 = this.resampleFunctionalConnectivityAll(ipr.sl_fc, sall);
            
            % expensive to build, so save  
            if ipr.save
                save([this.fileprefix('n', ipr.n) '.mat'], 'fc1')
            end
            
            %surf(fc1)
            %hold on
            h = figure;
            if ipr.z
                fc1 = atanh(fc1);
            end
            imagesc(fc1)
            set(gca, 'FontSize', 32)
            if isempty(ipr.colormap) || strcmpi(ipr.colormap, 'fccolormap')
                ld = load(fullfile(this.registry.workpath, 'FCColormap.mat'), 'FCColormap');
                colormap(gca, ld.FCColormap)
            else
                colormap(ipr.colormap)
            end
            if ipr.flip_colormap
                colormap(h, flipud(colormap(h)))
            end
            if ~isempty(ipr.caxis) 
                caxis(ipr.caxis)
            end
            colorbar
            %hold off
            h1 = gca; 
            h1.XAxis.TickLength = [0 0];
            h1.YAxis.TickLength = [0 0];
            set(gca, 'XTickLabel', [])
            set(gca, 'YTickLabel', [])
            %xlabel('voxels $\mathbf{x}$', 'Interpreter', 'latex', 'FontSize', 50)
            %ylabel('voxels $\mathbf{x}$', 'Interpreter', 'latex', 'FontSize', 50)  
            %this.savefig_highres(h);
            if ipr.save
                this.savefig(h, 'fileprefix', this.fileprefix('n', ipr.n));
            end
        end
        function fc1 = visualize_free_energy(this, varargin)
            %% visualizes Gramian for Carl's gray-matter mask.
            %
            % uhit coord	RSN (15)          ortho   ortho
            % 
            % 6:27		    smd                         1:22 
            % 28:31		    smi                        23:26
            % 32:50		    con                        27:45
            % 51:61		    aud                        46:56
            % 62:101		dmn                        57:96
            % 102:104		pmn                        97:99
            % 105:124		vis                       100:119
            % 125:139		fpn                       120:134
            % 140:147		sal                       135:142
            % 148:152		van                       143:147
            % 180:196		dan *             153:169 148:164
            % 153:167		bga               170:184 165:179
            % 168:179		tha               185:196 180:191
            % 197:198		amygdala    } MTL 197:201 192:196
            % 199:201		hippocampus }
            % 1:5			reward *                  197:201
            % * indicates relocation

            ip = inputParser;
            addOptional(ip, 'free_energy', [], @isnumeric)
            addParameter(ip, 'fileprefix', '', @ischar)
            addParameter(ip, 'caxis', [], @isnumeric)
            addParameter(ip, 'colormap', 'jet', @ischar)
            addParameter(ip, 'flip_colormap', false, @islogical)
            addParameter(ip, 'save', true, @islogical)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            bb300 = mlpetersen.BigBrain300('711-2B');
            %s300 = bb300.imagingContextSampling(301);
            sall = bb300.imagingContextSampling(20*300);
            fc1 = this.resampleFunctionalConnectivityReals(ipr.free_energy, sall);
            
            % expensive to build, so save  
            if ipr.save
                save(sprintf('%s_visualize_free_energy%s.mat', this.patientid, this.registry.tag), 'fc1')
            end
            
            %surf(fc1)
            %hold on
            h = figure;
            imagesc(fc1)
            set(gca, 'FontSize', 32)
            if isempty(ipr.colormap) || strcmpi(ipr.colormap, 'fccolormap')
                ld = load(fullfile(this.registry.workpath, 'FCColormap.mat'), 'FCColormap');
                colormap(gca, ld.FCColormap)
            else
                colormap(ipr.colormap)
            end
            if ipr.flip_colormap
                colormap(h, flipud(colormap(h)))
            end
            if ~isempty(ipr.caxis) 
                caxis(ipr.caxis)
            end
            colorbar
            %hold off
            h1 = gca; 
            h1.XAxis.TickLength = [0 0];
            h1.YAxis.TickLength = [0 0];
            set(gca, 'XTickLabel', [])
            set(gca, 'YTickLabel', [])
            %xlabel('voxels $\mathbf{x}$', 'Interpreter', 'latex', 'FontSize', 50)
            %ylabel('voxels $\mathbf{x}$', 'Interpreter', 'latex', 'FontSize', 50)  
            %this.savefig_highres(h);
            if ipr.save
                this.savefig(h, 'fileprefix', ipr.fileprefix);
            end
        end
        function pi_energy = visualize_pi_energy(this, varargin)
            ip = inputParser;
            addRequired(ip, 'fc1', @isnumeric)
            addParameter(ip, 'save', true, @islogical)
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
            if ipr.save
                this.savefig(h);
            end
        end
        
 		function this = SymmetricAFC(varargin)
 			%% SymmetricAFC
            %  @GMonly restricts to Carl's GM mask.
            %  @betaT is inverse temperature, default := 1.

 			this = this@mlafc.AFC(varargin{:});            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'patientid', '', @ischar)
            addParameter(ip, 'patientdir', pwd, @isfolder)
            addParameter(ip, 'similarityTag', '', @ischar)
            addParameter(ip, 'GMonly', true, @islogical)
            addParameter(ip, 'betaT', 1, @isscalar)
            addParameter(ip, 'tag', '_GMonly', @ischar)
            addParameter(ip, 'productTag', '_betaT1', @ischar)
            addParameter(ip, 'Nframes', [], @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.registry_.min_num_vox = 1;
            this.registry_.tanh_sandwich = true;
            this.registry_.similarityTag = ipr.similarityTag;
            this.registry_.tag = ipr.tag;
            this.registry_.productTag = ipr.productTag;
            this.registry_.GMonly = ipr.GMonly;
            this.registry.betaT = ipr.betaT;
            this.patientid_ = ipr.patientid;
            this.patientdir = ipr.patientdir;
            this.Nframes_ = ipr.Nframes;
            
            if isempty(this.sl_fc_mean_)
                this = load(this);
            end
        end        
    end 
    
    %% PROTECTED
    
    properties (Access = protected)
        Nframes_
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
    end
    
    %% HIDDEN
  
    methods (Static, Hidden)
    end
    
    methods (Hidden)
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

