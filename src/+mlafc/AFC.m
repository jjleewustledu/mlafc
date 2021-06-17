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
        patientid
 		registry
    end
    
    properties
        afc_map
        patientdir
    end
    
    methods (Static)
        function abs_Delta_fc = abs_Delta(fc, fc0)
            fcz = mlafc.AFC.atanh(double(fc));
            fcz0 = mlafc.AFC.atanh(double(fc0));
            abs_Delta_fc = tanh(abs(fcz - fcz0));
        end
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
        function ic = applyxfm(varargin)
            
            icref_ = fullfile(getenv('REFDIR'), '711-2B_111.nii.gz');
                
            ip = inputParser;
            addRequired(ip, 'ic0')
            addRequired(ip, 'mat', @isfile)
            addParameter(ip, 'out_fqfp', '', @ischar)
            addParameter(ip, 'icref', icref_)
            addParameter(ip, 'interp', 'trilinear', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            ipr.ic0 = mlfourd.ImagingContext2(ipr.ic0);
            ipr.icref = mlfourd.ImagingContext2(ipr.icref);
            
            bin = fullfile(getenv('FSLDIR'), 'bin', 'flirt');
            mlbash(sprintf('%s -in %s -applyxfm -init %s -out %s -paddingsize 0.0 -interp %s -ref %s', ...
                bin, ipr.ic0.fqfp, ipr.mat, ipr.out_fqfp, ipr.interp, ipr.icref.fqfp))
            %flirt -in PT15_4_seg_111_nopriorsurg.nii.gz -applyxfm -init PT15_4_FLIRT_on_mpr1_111_brain.mat -out PT15_4_seg_111_nopriorsurg_on_mpr_111.nii.gz -paddingsize 0.0 -interp nearestneighbour -ref PT15_mpr1_on_TRIO_Y_NDC_111.nii.gz
            ic = mlfourd.ImagingContext2([ipr.out_fqfp '.nii.gz']);
            ic.selectImagingFormatTool()
        end
        function x = atanh(x)
            %% avoids singularities for x -> {-1, 1}
            
            x = min(x, 1 - eps);
            x = max(x, -(1 - eps));
            x = atanh(x);
        end
        function ic = bet(varargin)
            ip = inputParser;
            addRequired(ip, 'obj', @(x) ~isempty(x))
            addOptional(ip, 'betFrac', 0.5, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            ic0 = mlafc.AFC.imagingContextNifti(ipr.obj);
            if ~isfile(ic0.fqfilename)
                ic0.save()
            end
            
            bin = fullfile(getenv('FSLDIR'), 'bin', 'bet');
            if contains(ic0.fileprefix, 'mpr')
                t2_fqfp = strrep(ic0.fqfp, 'mpr1', 't2w');
                assert(isfile([t2_fqfp '.nii.gz']))
                mlbash(sprintf('%s %s %s_brain -A2 %s -f %g -g 0 -m', bin, ic0.fqfp, ic0.fqfp, t2_fqfp, ipr.betFrac))
                BET = fullfile(ic0.filepath, 'BET', '');
                ensuredir(BET)
                movefile(fullfile(ic0.filepath, '*_mask.*'), BET, 'f')
                movefile(fullfile(ic0.filepath, '*_mesh.*'), BET, 'f')
            else
                mlbash(sprintf('%s %s %s_brain -R -f %i -g 0 -m', bin, ic0.fqfp, ic0.fqfp, ipr.betFrac))
            end
            ic = mlfourd.ImagingContext2([ic0.fqfp '_brain.nii.gz']);
            ic.selectImagingFormatTool()
        end
        function ic = betZ(varargin)
            ip = inputParser;
            addRequired(ip, 'obj', @(x) ~isempty(x))
            addOptional(ip, 'betFrac', 0.5, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            ic0 = mlafc.AFC.imagingContextNifti(ipr.obj);
            if ~isfile(ic0.fqfilename)
                ic0.save()
            end
            
            bin = fullfile(getenv('FSLDIR'), 'bin', 'bet');
            mlbash(sprintf('%s %s %s_brain -Z -f %i -g 0 -m', bin, ic0.fqfp, ic0.fqfp, ipr.betFrac))
            ic = mlfourd.ImagingContext2([ic0.fqfp '_brain.nii.gz']);
            ic.selectImagingFormatTool()
        end
        function calc_jsdiv(varargin)
            ip = inputParser;
            addOptional(ip, 'toglob', 'PT*', @ischar)
            addParameter(ip, 'outTag', '_softmax', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            for g = globFoldersT(ipr.toglob)
                pwd0 = pushd(g{1});
                
                segmentation = globT([g{1} '_*_segmentation_final_111.nii.gz']);
                seg = mlfourd.ImagingContext2(segmentation{1});
                afc_prob = mlfourd.ImagingContext2([g{1} ipr.outTag '_111.nii.gz']);
                gm = mlfourd.ImagingContext2(fullfile(getenv('REFDIR'), 'gm3d_111.nii.gz')); % no cerebellum
%                gm = gm.masked(double(afc_prob.numgt(0.008))); % 0.008 is the left tail of histograms                
                jsdiv = afc_prob.jsdiv(seg, gm);
                fprintf('jsdiv(%s) = %g\n',afc_prob.fileprefix, jsdiv)
                
                popd(pwd0)
            end
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
            addParameter(ip, 'ylabel', 'probability configuration')
            addParameter(ip, 'interpreter', 'tex', @ischar)
            addParameter(ip, 'transform', [])
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
                labels = {'2 (Ia)' '5 (IV)' '1 (Ia)' '6 (IV)' '4 (IIa)' '7 (IV)' '3 (Ia)'};  
                ordering = [3 1 7 5 2 4 6]; % pt28, pt15, pt36, pt34, pt26, pt29, p35
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
        function Delta_fc = Delta(fc, fc0)
            Delta_fc = tanh(mlafc.AFC.atanh(double(fc)) - mlafc.AFC.atanh(double(fc0)));
        end
        function img = flip12(img)
            img = flip(flip(img, 1), 2);
        end
        function [ic,matfn] = flirt(dof, ic0, icref, out_fqfp)
            assert(isnumeric(dof))
            bin = fullfile(getenv('FSLDIR'), 'bin', 'flirt');
            mlbash(sprintf('%s -in %s -ref %s -out %s -omat %s.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof %i -interp trilinear', bin, ic0.fqfp, icref.fqfp, out_fqfp, out_fqfp, dof))
            %flirt -in PT15_4_FLIRT_111_brain.nii.gz -ref PT15_mpr1_111_brain.nii.gz -out PT15_4_FLIRT_on_mpr1_111_brain.nii.gz -omat PT15_4_FLIRT_on_mpr1_111_brain.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear
            ic = mlfourd.ImagingContext2([out_fqfp '.nii.gz']);
            ic.selectImagingFormatTool()
            matfn = [out_fqfp '.mat'];
        end
        function flirt_to_MNI152(varargin)
            %  @param required nii is a NIfTI filename.
            %  @param interp is in {trilinear,nearestneighbour,sinc,spline}.
            
            ip = inputParser;
            addRequired(ip, 'nii', @isfile)
            addParameter(ip, 'fpout', '', @ischar)
            addParameter(ip, 'interp', 'trilinear', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            xfm = fullfile(getenv('REFDIR'), '711-2B_333_on_MNI152_T1_1mm.mat');
            if isempty(ipr.fpout)
                fnout = [myfileprefix(ipr.nii) '_on_MNI152.nii.gz'];
            else
                fnout = [ipr.fpout '.nii.gz'];
            end
            ref = fullfile(getenv('REFDIR'), '711-2B_333_on_MNI152_T1_1mm.nii.gz');
            mlbash(sprintf( ...
                '/usr/local/fsl/bin/flirt -in %s -applyxfm -init %s -out %s -paddingsize 0.0 -interp %s -ref %s', ...
                ipr.nii, xfm, fnout, ipr.interp, ref))
        end
        function ic = imagingContext_mask()
            reg = mlafc.AFCRegistry.instance();
            pth = fullfile(getenv('MLPDIR'), 'Reference_Images', '');
            if reg.GMonly
                ic = mlfourd.ImagingContext2( ...
                    fullfile(pth, 'N21_aparc+aseg_GMctx_on_711-2V_333_avg_zlt0.5_gAAmask_v1.4dfp.hdr'));
            else                
                ic = mlfourd.ImagingContext2( ...
                    fullfile(pth, 'glm_atlas_mask_333.4dfp.hdr'));
            end
        end
        function ic = imagingContextNifti(varargin)
            ip = inputParser;
            addRequired(ip, 'obj', @(x) ~isempty(x))
            parse(ip, varargin{:})
            ipr = ip.Results;
            ipr.obj = mlfourd.ImagingContext2(ipr.obj);
            ipr.obj = mlfourd.ImagingContext2(ipr.obj.nifti);
            fp = ipr.obj.fileprefix;
            
            if contains(fp, '111')
                nii = mlfourd.ImagingFormatContext(fullfile(getenv('REFDIR'), '711-2B_111.nii.gz'));
            elseif contains(fp, '222')
                nii = mlfourd.ImagingFormatContext(fullfile(getenv('REFDIR'), '711-2B_222.nii.gz'));
            elseif contains(fp, '333')
                nii = mlfourd.ImagingFormatContext(fullfile(getenv('REFDIR'), '711-2B_333.nii.gz'));
            else
                ic = ipr.obj;
                return
            end
            nii.img = ipr.obj.nifti.img;
            nii.filepath = ipr.obj.filepath;
            nii.fileprefix = ipr.obj.fileprefix;
            ic = mlfourd.ImagingContext2(nii);
        end
        function d = kldiv(corrP, corrQ)
            %% KLDIV estimates divergence between the prob of FC for data and the prob of FC for a model
            %  using complete asymmetric correlation matrices without downsampling.
            %  @param corrP describes spheres-fiber x base for data.
            %  @param corrQ describes spheres-fiber x base for model.
            %  @return KL divergence in vector ~ 1 x base in nats.
            
            P = exp(corrP); % prob of corr matrix for data
            P(isnan(P)) = eps;
            P(P < eps) = eps; % manage round-off
            for ib = 1:size(P,2)
                P(:,ib) = P(:,ib) / sum(P(:,ib),1); % ~ prob normalized to each base voxel
            end
            
            Q = exp(corrQ); % prob of corr matrix for model
            Q(isnan(Q)) = eps;
            Q(Q < eps) = eps; % manage round-off
            for ib = 1:size(Q,2)
                Q(:,ib) = Q(:,ib) / sum(Q(:,ib),1); % ~ prob normalized to each base voxel
            end
            
            d = sum(P .* (log(P) - log(Q)), 1); % marginalize fibers
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
            obj = mlafc.SymmetricAFC.maskArrToFullArr(obj);
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
            ifc = mlfourd.ImagingFormatContext(fullfile(getenv('REFDIR'), '711-2B_333.nii.gz'));
            ifc.img = img;
            ifc.filepath = pwd;
            ifc.fileprefix = myfileprefix(matname);
            ic = mlfourd.ImagingContext2(ifc);
        end    
        function kl_smpl = probToKLSample(q)
            %  @param q is numeric, a probability for each element.
            %  @return kl is numeric, the integrand of K-L divergence from the uniform distrib. for each element.
            %          uniform distrib. ~ exp(-| Delta() |) ~ 1 / N, with Delta := 0.
            
            reg = mlafc.AFCRegistry.instance();
            N = reg.ref_count;
            
            % KL := -E_p log(q/p), p uniform ~ 1/500, q is approximating.
            kl_smpl = -log(N*q); % just one sample for q
        end
        function vec = imgread(fn)
            assert(isfile(fn))
            [pth,fp,x] = myfileparts(fn);
            if lstrfind(x, '4dfp')
                vec = mlperceptron.Fourdfp.Read4dfp(fullfile(pth, [fp '.4dfp.img']));
                return
            end
            if lstrfind(x, '.nii')
                system(sprintf('niftigz_4dfp -4 %s%s %s.4dfp.hdr', fullfile(pth, fp), x, fullfile(pth, fp)))
                vec = mlperceptron.Fourdfp.Read4dfp(fullfile(pth, [fp '.4dfp.img']));
                deleteExisting(fullfile(pth, [fp '.4dfp.*']))
            end
        end
    end
    
    methods
        
        %% GET
        
        function g = get.patientid(this)
            if ~isempty(this.patientid_)
                g = this.patientid_;
                return
            end
            pid = mybasename(this.patientdir);
            if contains(lower(pid), 'pt')
                g = pid;
            else
                g = '';
            end
        end
        function this = set.patientid(this, s)
            assert(ischar(s))
            this.patientid_ = s;
        end
        function g = get.registry(this)
            g = this.registry_;
        end
        
        %%
		  
 		function this = AFC(varargin)
            import mlperceptron.*;
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'sphere_radius', 1, @isnumeric)
            addParameter(ip, 'grid_spacing', 1, @isnumeric)
            addParameter(ip, 'min_num_vox', 1, @isnumeric)
            addParameter(ip, 'tanh_sandwich', true, @islogical)
            addParameter(ip, 'ref_count', 500, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.registry_ = mlafc.AFCRegistry.instance('initialize'); 
            this.registry_.sphere_radius = ipr.sphere_radius;
            this.registry_.grid_spacing = ipr.grid_spacing;
            this.registry_.min_num_vox = ipr.min_num_vox;
            this.registry_.tanh_sandwich = ipr.tanh_sandwich;
            this.registry_.ref_count = ipr.ref_count;
            this.registry_.tag = '';
            
            this.glmatl_ = PerceptronRegistry.read_glm_atlas_mask();
            this.glmatl_(find(this.glmatl_)) = 1;  
            this.glmmsk_indices_ = find(this.glmatl_);
            GMctx = PerceptronRegistry.read_N21_aparc_aseg_GMctx(); % 147456 x 1
            this.GMmsk_ = GMctx;
            this.GMmsk_(find(GMctx)) = 1;
            this.GMmsk_for_glm_ = GMctx(find(this.glmatl_)); %#ok<*FNDSB>             
 		end
    end 
    
    %% PROTECTED
    
    properties % (Access = protected)
        glmatl_         % 147456 x 1 ~ 48*64*48 x 1
        glmmsk_indices_ % 65549  x 1 ~ |parenchyma|
        GMmsk_          % 147456 x 1 ~ 48*64*48 x 1
        GMmsk_for_glm_  % 65549  x 1 ~ |parenchyma|
        patientid_
        registry_
        sl_fc_mean_     % 2825 x 18611 ~ |sv_| x |grey matter|
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

