classdef AFCRegistry < handle & mlpark.ParkRegistry
	%% AFCREGISTRY  

	%  $Revision$
 	%  was created 31-May-2019 16:35:57 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlafc/src/+mlafc.
 	%% It was developed on Matlab 9.5.0.1067069 (R2018b) Update 4 for MACI64.  Copyright 2019 John Joowon Lee.
 	
	properties (Dependent)
        Hacker_Data_ALL
        MLPAFC_dir
        MLP_GTM_100
        mlp_rmse_con100
        perceptron_resid_mat
        ref_count
        ref_resid_mat
 		sl_fc_gsp_mat
        sl_fc_mean_mat
 	end

    methods (Static)
        function [nx,ny,nz,n3d] = atlas_dims(varargin)
            [nx,ny,nz,n3d] = mlperceptron.PerceptronRegistry.atlas_dims(varargin{:});
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
        
        function g = get.Hacker_Data_ALL(this)
            g = fullfile(this.root, 'data', 'nil-bluearc', 'corbetta', 'Hacker', 'Data', 'ALL', '');
        end
        function g = get.MLPAFC_dir(this)
            g = fullfile(this.parkhome, 'MLPAFC', '');
        end
        function g = get.MLP_GTM_100(this)
            if isempty(this.MLP_GTM_100_)                
                load(fullfile(this.MLPAFC_dir, 'MLP_GTM_100.mat'), 'MLP_GTM_100')
                this.MLP_GTM_100_ = MLP_GTM_100;  %#ok<PROP>
                clear('MLP_GTM_100')
            end
            g = this.MLP_GTM_100_;
        end
        function g = get.mlp_rmse_con100(this)
            if isempty(this.mlp_rmse_con100_)                
                load(fullfile(this.MLPAFC_dir, 'mlp_rmse_con100.mat'), 'mlp_rmse_con100')
                this.mlp_rmse_con100_ = mlp_rmse_con100;  %#ok<PROP>
                clear('mlp_rmse_con100')
            end
            g = this.mlp_rmse_con100_;
        end
        function g = get.perceptron_resid_mat(~)
            %g = '_faln_dbnd_xr3d_uwrp_atl_uout_resid.mat';
            g = '_faln_dbnd_xr3d_atl_g7_bpss_resid.mat';
        end
        function g = get.ref_count(~)
            g = 100;
        end
        function g = get.ref_resid_mat(~)
            g = '_faln_dbnd_xr3d_atl_g7_bpss_resid.mat';
        end
        function g = get.sl_fc_gsp_mat(this)
            g = fullfile(getenv('WORK'), 'sl_fc_gsp.mat');
        end
        function g = get.sl_fc_mean_mat(this)
            g = fullfile(getenv('WORK'), 'sl_fc_mean.mat');
        end
    end
    
    %% PROTECTED
    
    properties (Access = protected)
        MLP_GTM_100_
        mlp_rmse_con100_
    end
    
	methods (Access = protected)
		  
 		function this = AFCRegistry(varargin)
 			this = this@mlpark.ParkRegistry(varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

