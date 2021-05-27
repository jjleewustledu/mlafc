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
        sv
    end
    
    properties
        afc_map
        patientdir
    end
    
    methods (Static)
        function ic = applyxfm(ic0, mat, out_fqfp, icref)
            bin = fullfile(getenv('FSLDIR'), 'bin', 'flirt');
            mlbash(sprintf('%s -in %s -applyxfm -init %s -out %s -paddingsize 0.0 -interp nearestneighbour -ref %s', bin, ic0.fqfp, mat, out_fqfp, icref.fqfp))
            %flirt -in PT15_4_seg_111_nopriorsurg.nii.gz -applyxfm -init PT15_4_FLIRT_on_mpr1_111_brain.mat -out PT15_4_seg_111_nopriorsurg_on_mpr_111.nii.gz -paddingsize 0.0 -interp nearestneighbour -ref PT15_mpr1_on_TRIO_Y_NDC_111.nii.gz
            ic = mlfourd.ImagingContext2([out_fqfp '.nii.gz']);
            ic.selectImagingFormatTool()
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
        function g = get.sv(this)
            g = this.sv_;
        end
        
        %%
		  
 		function this = AFC(varargin)
            import mlperceptron.*;
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'sphere_radius', 2, @isnumeric)
            addParameter(ip, 'grid_spacing', 3, @isnumeric)
            addParameter(ip, 'min_num_vox', 1, @isnumeric)
            addParameter(ip, 'tanh_sandwich', false, @islogical)
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
            this.GMmsk_for_glm_ = GMctx(find(this.glmatl_)); %#ok<*FNDSB>             
 		end
    end 
    
    %% PROTECTED
    
    properties % (Access = protected)
        glmatl_         % 147456 x 1 ~ 48*64*48 x 1
        glmmsk_indices_ % 65549  x 1 ~ |parenchyma|
        GMmsk_for_glm_  % 65549  x 1 ~ |parenchyma|
        patientid_
        registry_
        sl_fc_mean_     % 2825 x 18611 ~ |sv_| x |grey matter|
        sv_             % 1    x 2825  cells with variable size [(x > 1) 1]
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

