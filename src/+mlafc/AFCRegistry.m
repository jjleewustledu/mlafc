classdef AFCRegistry < handle
	%% AFCREGISTRY  

	%  $Revision$
 	%  was created 31-May-2019 16:35:57 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlafc/src/+mlafc.
 	%% It was developed on Matlab 9.5.0.1067069 (R2018b) Update 4 for MACI64.  Copyright 2019 John Joowon Lee.
 	
    properties
        atlasCode = 333
        betaT = 1
        bold_suffix = '_faln_dbnd_xr3d_atl_g7_bpss_resid'
        GMonly = true
        grid_spacing = 1
        min_num_vox = 1
        productTag = '_betaT1'
        ref_count = 500
        similarityTag = ''
        sphere_radius = 1
        tag = '_GMonly'
        tanh_sandwich = true
    end
    
	properties (Dependent)
        parkhome
        gtm500_dir
        gtm500_ids
        
        afc_map_mat
        fileTag
        Hacker_Data_ALL
        N_BOLD
        perceptron_uout_resid_mat
        perceptron_resid_mat
        ref_resid_mat
        sl_fc_gsp_mat
 		sl_fc_gsp_sum_prob_mat
        sl_fc_mean_mat
        sl_fc_std_mat
        tanh_tag
        workpath
        Z_E_mat
        Z_G_mat
 	end

    methods (Static)
        function [nx,ny,nz,n3d] = atlas_dims(varargin)
            [nx,ny,nz,n3d] = mlperceptron.PerceptronRegistry.atlas_dims(varargin{:});
        end
        function n3d            = atlas_numel(varargin)
            n3d = mlperceptron.PerceptronRegistry.atlas_numel(varargin{:});
        end
        function ns             = atlas_dims_vec(varargin)
            ns = mlperceptron.PerceptronRegistry.atlas_dims_vec(varargin{:});
        end
        function ic             = segmentation(pid)
            assert(ischar(pid))
            pid = upper(pid);
            pwd0 = pushd(fullfile('/data/nil-bluearc/shimony/jjlee/FocalEpilepsy/Emily/Segmentations_06_and_07_2020', pid));
            fn = glob([pid '*eg*_111_to_333.nii.gz']);
            ic = mlfourd.ImagingContext2(fn);
            popd(pwd0)
        end
        
        function this = instance(qualifier)
            %% INSTANCE uses string qualifiers to implement registry behavior that
            %  requires access to the persistent uniqueInstance
            %  @optional qualifier := 'initialize' | ''.

            persistent uniqueInstance
            
            if (exist('qualifier','var') && ischar(qualifier))
                if (strcmp(qualifier, 'initialize'))
                    uniqueInstance = [];
                end
            end
            
            if (isempty(uniqueInstance))
                this = mlafc.AFCRegistry();
                uniqueInstance = this;
            else
                this = uniqueInstance;
            end
        end
    end 
    
    methods 
        
        %% GET
        
        function g = get.afc_map_mat(~)
            g = ['afc_map_' datestr(now, 30)];
        end
        function g = get.fileTag(this)
            g = sprintf('_radius%i_stride%i%s%s', this.sphere_radius, this.grid_spacing, this.similarityTag, this.tag);
        end
        function g = get.gtm500_dir(this)
            g = fullfile(this.parkhome, 'GTM500Perceptron', ''); 
        end
        function g = get.gtm500_ids(this)
            if isempty(this.gtm500_ids_)
                this.gtm500_ids_ = importdata(fullfile(this.gtm500_dir, 'GTM500.lst'));
            end
            g = this.gtm500_ids_;
        end
        function g = get.Hacker_Data_ALL(~)
            g = fullfile('/', 'data', 'nil-bluearc', 'corbetta', 'Hacker', 'Data', 'ALL', '');
        end        
        function g = get.N_BOLD(this)
            if this.GMonly
                g = 18611;
            else
                g = 65549;
            end
        end
        function g = get.parkhome(~)
            g = '/data/nil-bluearc/shimony/Park';
        end
        function g = get.perceptron_uout_resid_mat(~)
            g = '_faln_dbnd_xr3d_uwrp_atl_uout_resid.mat';
        end
        function g = get.perceptron_resid_mat(~)
            g = '_faln_dbnd_xr3d_atl_g7_bpss_resid.mat';
        end
        function g = get.ref_resid_mat(~)
            g = '_faln_dbnd_xr3d_atl_g7_bpss_resid.mat';
        end
        function g = get.sl_fc_gsp_mat(this)
            g = fullfile(this.workpath, ...
                sprintf('sl_fc_gsp_radius%i_stride%i_N%i%s%s.mat', ...
                this.sphere_radius, this.grid_spacing, this.ref_count, this.tanh_tag, this.tag));
        end
        function g = get.sl_fc_gsp_sum_prob_mat(this)
            g = fullfile(this.workpath, ...
                sprintf('sl_fc_gsp_sum_prob_radius%i_stride%i_N%i%s%s%s.mat', ...
                this.sphere_radius, this.grid_spacing, this.ref_count, this.tanh_tag, this.tag, this.similarityTag));
        end
        function g = get.sl_fc_mean_mat(this)
            g = fullfile(this.workpath, ...
                sprintf('sl_fc_mean_radius%i_stride%i_N%i%s%s.mat', ...
                this.sphere_radius, this.grid_spacing, this.ref_count, this.tanh_tag, this.tag));
        end
        function g = get.sl_fc_std_mat(this)
            g = fullfile(this.workpath, ...
                sprintf('sl_fc_std_radius%i_stride%i_N%i%s%s.mat', ...
                this.sphere_radius, this.grid_spacing, this.ref_count, this.tanh_tag, this.tag));
        end
        function g = get.Z_E_mat(this)
            g = fullfile(pwd, ...
                sprintf('Z_E_radius%i_stride%i_N%i%s.mat', ...
                this.sphere_radius, this.grid_spacing, this.ref_count, this.tag));
        end
        function g = get.Z_G_mat(this)
            g = fullfile(this.workpath, ...
                sprintf('Z_G_radius%i_stride%i_N%i%s.mat', ...
                this.sphere_radius, this.grid_spacing, this.ref_count, this.tag));
        end
        function g = get.tanh_tag(this)
            if this.tanh_sandwich
                g = '_tanh';
            else 
                g = '';
            end
        end
        function g = get.workpath(~)
            if contains(hostname, 'precuneal')
                g = fullfile(getenv('HOME'), 'Box', 'Leuthardt_Epilepsy_Project', 'Segmentations_06_and_07_2020', '');
            else
                g = '/data/nil-bluearc/shimony/jjlee/FocalEpilepsy';
            end
        end
        
        %%        
        
        function m = sl_fc_gsp_ref_mat(this, ref)
            assert(isscalar(ref))
            m = fullfile(this.workpath, ...
                sprintf('sl_fc_gsp%i_radius%i_stride%i_N%i%s.mat', ...
                ref, this.sphere_radius, this.grid_spacing, this.ref_count, this.tag));
        end
    end
    
    %% PROTECTED
    
    properties (Access = protected)
        gtm500_ids_
    end
    
	methods (Access = protected)		  
 		function this = AFCRegistry(varargin)
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

