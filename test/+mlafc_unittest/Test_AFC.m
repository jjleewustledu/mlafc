classdef Test_AFC < matlab.unittest.TestCase
	%% TEST_AFC 

	%  Usage:  >> results = run(mlafc_unittest.Test_AFC)
 	%          >> result  = run(mlafc_unittest.Test_AFC, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 29-May-2019 16:35:03 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlafc/test/+mlafc_unittest.
 	%% It was developed on Matlab 9.5.0.1067069 (R2018b) Update 4 for MACI64.  Copyright 2019 John Joowon Lee.
 	
	properties
 		registry
 		testObj
        WORK 
        %WORK = '/data/shimony/shimony4/jjlee/FocalEpilepsy'
    end
    
    properties (Dependent)
        sv_ref
    end
    
    methods
        function g = get.sv_ref(this)
            load(fullfile(this.WORK, 'SL_fMRI_initialization.mat'), 'sv')
            g = sv;
        end
    end

	methods (Test)
        function test_mapOfSpheres(~)
            for sr = [3 2 1]
                jafc = mlafc.SymmetricAFC( ...
                    'sphere_radius', sr, 'grid_spacing', sr);
                ic = jafc.mapOfSpheres();
                ic.fsleyes()
            end
        end
        function test_radius(this)
            setenv('WORK', this.WORK)
            cd(this.WORK)
            
            tic
            jafc = mlafc.SymmetricAFC( ...
                'sphere_radius', 1, ...
                'grid_spacing', 1, ...
                'ref_count', 10);
            jafc = jafc.explore_fc();
            disp(jafc)
            ic = mlafc.SymmetricAFC.slfcMatToIC('sl_fc_mean_radius1_stride1_N10_tanh_GMonly.mat', 'sl_fc_mean');
            toc
            ic.save()
            ic.fsleyes()
            
            tic
            jafc = mlafc.SymmetricAFC( ...
                'sphere_radius', 2, ...
                'grid_spacing', 3, ...
                'ref_count', 10);
            jafc = jafc.explore_fc();
            disp(jafc)
            ic = mlafc.SymmetricAFC.slfcMatToIC('sl_fc_mean_radius2_stride3_N10_tanh_GMonly.mat', 'sl_fc_mean');
            toc
            ic.save()
            ic.fsleyes()
        end
        function test_explore_fc(this)
            setenv('WORK', this.WORK)
            cd(this.WORK)
            jafc = mlafc.SymmetricAFC( ...
                'sphere_radius', 1, ...
                'grid_spacing', 1, ...
                'ref_count', 500);
            jafc = jafc.explore_fc();
            disp(jafc)
            ic = mlafc.SymmetricAFC.slfcMatToIC('sl_fc_mean_radius1_stride1_N500_tanh_GMonly.mat', 'sl_fc_mean');
            ic.save()
            ic.fsleyes()
        end
        function test_explore_fc2(this)
            setenv('WORK', this.WORK)
            cd(this.WORK)
            jafc = mlafc.SymmetricAFC( ...
                'sphere_radius', 1, ...
                'grid_spacing', 1, ...
                'ref_count', 500);
            jafc = jafc.explore_fc2();
            disp(jafc)
            ic = mlafc.SymmetricAFC.slfcMatToIC('sl_fc_std_radius1_stride1_N500_tanh_GMonly.mat', 'sl_fc_std');
            ic.save()
            ic.fsleyes()
        end
        function test_make_sl_fc_intermediates(this)
            setenv('WORK', this.WORK)
            cd(this.WORK)
            pts = globFoldersT('PT*');
            for ip = 1:length(pts)
                if strcmp(pts{ip}, 'PTmore')
                    continue
                end
                try
                    pdir = fullfile(getenv('WORK'), pts{ip}, '');
                    pwd0 = pushd(pdir);
                    jafc = mlafc.SymmetricAFC( ...
                        'sphere_radius', 2, ...
                        'grid_spacing', 3, ...
                        'ref_count', 500);
                    jafc = jafc.make_sl_fc_intermediates(pwd, pts{ip});
                    popd(pwd0)
                catch ME
                    handwarning(ME)
                end
            end            
        end
        function test_makeSoftmax(this)
            setenv('WORK', this.WORK)
            cd(this.WORK)
            pts = globFoldersT('PT*');
            for ip = 1:length(pts)
                if strcmp(pts{ip}, 'PTmore')
                    continue
                end
                try
                    pdir = fullfile(getenv('WORK'), pts{ip}, '');
                    pwd0 = pushd(pdir);
                    jafc = mlafc.JohnsAFC( ...
                        'sphere_radius', 1, ...
                        'grid_spacing', 1, ...
                        'ref_count', 500);
                    jafc = jafc.makeSoftmax(pwd, pts{ip});
                    popd(pwd0)
                catch ME
                    handwarning(ME)
                end
            end            
        end
        function test_makeSoftmax1(this)
            setenv('WORK', this.WORK)
            cd(this.WORK)
            pts = globFoldersT('PT*');
            parfor ip = 1:length(pts)
                try
                    pdir = fullfile(getenv('WORK'), pts{ip}, '');
                    pwd0 = pushd(pdir);
                    jafc = mlafc.SymmetricAFC( ...
                        'betaT', 1, ...
                        'productTag', '_betaT1', ...
                        'Nframes', 1200);
                    jafc = jafc.makeSoftmax(pwd, pts{ip});
                    disp(jafc)
                    popd(pwd0)
                catch ME
                    handwarning(ME)
                end
            end            
        end
        function test_makeSoftmaxGsp(this)
            setenv('WORK', this.WORK)
            cd(this.WORK)
            pts = globFoldersT('Sub*_Ses1');
            parfor ip = 1:length(pts)
                try
                    pdir = fullfile(getenv('WORK'), pts{ip}, '');
                    pwd0 = pushd(pdir);
                    jafc = mlafc.SymmetricAFC( ...
                        'productTag', '_gsp');
                    jafc = jafc.makeSoftmax(pwd, pts{ip});
                    disp(jafc)
                    popd(pwd0)
                catch ME
                    handwarning(ME)
                end
            end            
        end
        function test_makeSoftmaxMsc(this)
            setenv('WORK', this.WORK)
            cd(this.WORK)
            pts = globFoldersT('MSC*');
            for ip = 1:length(pts)
                try
                    pdir = fullfile(getenv('WORK'), pts{ip}, '');
                    pwd0 = pushd(pdir);
                    jafc = mlafc.SymmetricAFC( ...
                        'sphere_radius', 2, ...
                        'grid_spacing', 3, ...
                        'ref_count', 500);
                    jafc = jafc.makeSoftmax(pwd, pts{ip}, 'Nframes', 2307);
                    jafc.product.save() % .4dfp.*
                    jafc.product.nifti.save() % .nii.gz
                    popd(pwd0)
                catch ME
                    handwarning(ME)
                end
            end            
        end
        function test_visualizePT(this)
            setenv('WORK', this.WORK)
            cd(this.WORK)
            ld = load('Z_G_radius1_stride1_N500_GMonly.mat', 'Z_G');
            Z_G = ld.Z_G;
            pts = globFoldersT('PT*');
            outmap = containers.Map('KeyType', 'char', 'ValueType', 'any');
            for ip = 1:length(pts)
                if strcmp(pts{ip}, 'PTmore')
                    continue
                end
                try
                    pdir = fullfile(getenv('WORK'), pts{ip}, '');
                    pwd0 = pushd(pdir);
                    outmap(pts{ip}) = mlafc.SymmetricAFC.visualizePT(pts{ip}, Z_G, 'KLSample', true, 'save', false);
                    popd(pwd0)
                catch ME
                    handwarning(ME)
                end
            end 
            this.report_out(outmap)
        end
        function test_fullArrToMaskArr(this)
            atl = mlfourd.ImagingFormatContext('711-2B_333.nii.gz');
            farr = reshape(flip(flip(atl.img,1),2), [1 48*64*48]);
            garr = mlafc.SymmetricAFC.fullArrToMaskArr(farr);
            farr = mlafc.SymmetricAFC.maskArrToFullArr(garr);
            atl.img = reshape(farr, [48 64 48]);
            atl.img = flip(flip(atl.img, 1), 2);
            atl.fileprefix = 'test';
            atl.fsleyes
        end
        function test_x333_to_xMaskArr(this)
            if ~this.registry.GMonly
                this.assertEqual(mlafc.SymmetricAFC.x333_to_xMaskArr([23 48 21]), 25878)
                this.assertEqual(mlafc.SymmetricAFC.x333_to_xMaskArr([36 48 21]), 25865)
                this.assertEqual(mlafc.SymmetricAFC.x333_to_xMaskArr([11 14 35]), 55267)
            end
        end
        function test_buildFiberBundles_contrast(this)
            setenv('WORK', this.WORK)
            cd(this.WORK)
            sl_fc_mean_ic = mlafc.SymmetricAFC.slfcMatToIC(this.registry.sl_fc_mean_mat, 'sl_fc_mean');
            pts = globFoldersT('PT*');
            for ip = 1:length(pts)
                if strcmp(pts{ip}, 'PTmore')
                    continue
                end
                try
                    pdir = fullfile(getenv('WORK'), pts{ip}, '');
                    pwd0 = pushd(pdir);
                    mlafc.SymmetricAFC.buildFiberBundles_contrast( ...
                        [pts{ip} '_sl_fc.mat'], ...
                        'sl_fc_mean_ic', sl_fc_mean_ic)
                    popd(pwd0)
                catch ME
                    handwarning(ME)
                end
            end  
        end
        
        function test_make_SLfMRI_init(this)
            setenv('WORK', this.WORK)
            cd(this.WORK)
            this.testObj.SL_fMRI_initialization()
        end
        function test_makeSoftmax_ori(this)
            setenv('WORK', this.WORK)
            cd(this.WORK)
            pts = {'PT15' 'PT26' 'PT28' 'PT34' 'PT35' 'PT36'};
            for ip = 1:length(pts)
                pdir = fullfile(getenv('WORK'), pts{ip}, '');
                pwd0 = pushd(pdir);
                testObj = this.testObj.makeSoftmax(pdir, pts{ip});
                testObj.product.fsleyes()
                testObj.product.save()
                testObj.product.nifti.save()
                popd(pwd0)
            end            
        end
        function test_makeDice(this)
            setenv('WORK', this.WORK)
            cd(this.WORK)
            testObj = mlafc.EmilysAFC();
            pts = {'PT15' 'PT26' 'PT28' 'PT34' 'PT35' 'PT36'};
            sims = zeros(1, length(pts));
            summary = '';
            for s = 0:3
                for ip = 1:length(pts)
                    pdir = fullfile(getenv('WORK'), pts{ip}, '');
                    sdir = fullfile(getenv('WORK'), 'Emily', 'Segmentations_06_and_07_2020', pts{ip}, '');
                    pwd0 = pushd(pdir);
                    sims(ip) = testObj.makeDice( ...
                        pdir, pts{ip}, ...
                        'afc_filename', fullfile(pdir, 'AFC_and_resections', 'AFC_and_resections_obj.mat'), ...
                        'feature_filename', fullfile(sdir, [pts{ip} '_seg_toBoldDeformed.nii.gz']), ...
                        'load_afc', true, ...
                        'Nsigma', s);
                    popd(pwd0)
                end
                save(sprintf('test_makeDice_%isd.mat', s), 'sims')
                for ip = 1:length(pts)
                    summary = [summary sprintf('%s %isd similarity -> %g\n', pts{ip}, s, sims(ip))]; %#ok<AGROW>
                end
            end    
            fprintf(summary)
        end
        function test_msc(this)
            obj = mlafc.AFCFromMat;
            obj = obj.makeMscMap( ...
                fullfile(getenv('WORK'), 'MSC_subject6', ''), 'MSC_subject6', ...
                'afc_filename', 'AFC_this_msc_subject6.mat');
            disp(obj)
        end
        function test_hotspot_corr(this)
            pdir = fullfile(getenv('WORK'), 'PT36', '');
            pwd0 = pushd(pdir);
            hotspot = this.testObj.hotspot_corr( ...
                pdir, 'PT36', ...
                'afc_filename', fullfile(pdir, 'AFC_and_resections', 'AFC_and_resections.mat'), ...
                'feature_filename', 'epileptiformECoG_mark333.4dfp.hdr', ...
                'load_afc', true); % afc.4dfp.hdr
            hotspot.fsleyes
            hotspot.save
            popd(pwd0)
        end
        function test_resection_corr(this)
            pdir = fullfile(getenv('WORK'), 'PT36', '');
            pwd0 = pushd(pdir);
            hotspot = this.testObj.resection_corr( ...
                pdir, 'PT36', ...
                'afc_filename', fullfile(pdir, 'AFC_and_resections', 'AFC_and_resections.mat'), ...
                'feature_filename', fullfile(pdir, 'PT36_Segmentation_333.nii'), ...
                'load_afc', true);
            hotspot.fsleyes
            hotspot.save
            popd(pwd0)
        end
        function test_ctor(this)
            disp(this.testObj)
        end
		function test_SL_fMRI_initialization(this)
 			import mlafc.*;
            this.testObj = this.testObj.SL_fMRI_initialization( ...
                'make_sl_fc_mean', false, 'make_sl_fc_gsp', false);
            sv = this.sv_ref;
            this.verifyEqual(size(this.testObj.sv_), size(sv))
            for s = length(sv)
                this.verifyEqual(this.testObj.sv_{s}, sv{s})
            end
        end
        function test_SL_AFC(this)
            % cp from Carl's repo to fullfile(getenv('WORK'), 'PT36', ''); du -hs ~ 66G
            this.testObj.SL_AFC( ...
                fullfile(getenv('WORK'), 'PT36', ''), 'PT36')
        end
        function test_figures(this)            
            this.testObj = this.testObj.makeDifferenceMap( ...
                fullfile(getenv('WORK'), 'PT36', ''), 'PT36', 'afc_filename', 'AFC_this_1400.mat');
            
            for n = 300:300:1200
                try
                    this.testObj.makeDifferenceMap( ...
                        fullfile(getenv('WORK'), 'PT36', ''), 'PT36', 'Nframes', n,  'afc_filename', sprintf('AFC_this_%i.mat', n));
                catch ME
                    handwarning(ME)
                end
            end
        end
        function test_meeting_20200117(this)
            pts = {'PT36' 'PT26' 'PT28' 'PT3' 'PT16' 'PT5' 'PT15' 'PT34' 'PT35'};
            for p = pts
                pdir = fullfile(getenv('WORK'), p{1}, '');
                pwd0 = pushd(pdir);
                this.testObj.makeROC( ...
                    pdir, p{1}, ...
                    'afc_filename', fullfile(pdir, 'AFC_and_resections', 'AFC_and_resections.mat'), ...
                    'feature_filename', fullfile(pdir, [p{1} '_Segmentation_333.nii']), ...
                    'load_afc', true);
                popd(pwd0)
            end
        end
        function test_meeting_20200501(this)
            pts = {'PT36' 'PT26' 'PT28' 'PT3' 'PT16' 'PT5' 'PT15' 'PT34' 'PT35'};
            sims = zeros(1, length(pts));
            for ip = 1:length(pts)
                pdir = fullfile(getenv('WORK'), pts{ip}, '');
                pwd0 = pushd(pdir);
                sims(ip) = this.testObj.makeDice( ...
                    pdir, pts{ip}, ...
                    'afc_filename', fullfile(pdir, 'AFC_and_resections', 'AFC_and_resections.mat'), ...
                    'feature_filename', fullfile(pdir, [pts{ip} '_Segmentation_333.nii']), ...
                    'load_afc', true);
                popd(pwd0)
            end
            for ip = 1:length(pts)
                fprintf('%s similarity -> %g\n', pts{ip}, sims(ip))
            end
        end
        function test_New_Epilepsy(this)
            newepil = '/data/nil-bluearc/shimony/jjlee/FocalEpilepsy';
            pwd0 = pushd(newepil);
            dt = mlsystem.DirTool('PT*');
            for i = 1:length(dt.fqdns)
                this.testObj.SL_AFC(dt.fqdns{i}, dt.dns{i})
            end
            this.testObj.SL_AFC(fullfile(newepil, '09_0009'), '09_0009')
            popd(pwd0)
        end
        function test_PT34(this)
            pdir = fullfile(getenv('WORK'), 'PT34', '');
            pwd0 = pushd(pdir);
            this.testObj.makeROC( ...
                pdir, 'PT34', ...
                'afc_filename', fullfile(pdir, 'AFC_and_resections', 'AFC_and_resections.mat'), ...
                'feature_filename', fullfile(pdir, 'PT34_Segmentation_333.nii'), ...
                'load_afc', true);
            popd(pwd0)
        end
        function test_PT3(this)
            pdir = fullfile(getenv('WORK'), 'PT3', '');
            pwd0 = pushd(pdir);
            this.testObj.makeROC( ...
                pdir, 'PT3', ...
                'afc_filename', fullfile(pdir, 'AFC_and_resections', 'AFC_and_resections.mat'), ...
                'feature_filename', fullfile(pdir, 'PT3_Segmentation_333.nii'), ...
                'load_afc', true);
            popd(pwd0)
        end
        function test_09_0009(this)
            newepil = '/data/nil-bluearc/shimony/jjlee/FocalEpilepsy';
            pwd0 = pushd(newepil);
            this.testObj.SL_AFC(fullfile(newepil, '09_0009'), '09_0009')
            popd(pwd0)            
        end
        function test_TLE(this)
            ltle = '/data/nil-bluearc/shimony/jjlee/TLE/LTLE';
            rtle = '/data/nil-bluearc/shimony/jjlee/TLE/RTLE';
            pwd0 = pushd(ltle);
            dt = mlsystem.DirTool('09_0*');
            for i = 1:length(dt.fqdns)
                this.testObj.SL_AFC(dt.fqdns{i}, dt.dns{i})
            end
            popd(pwd0)
            pwd0 = pushd(rtle);
            dt = mlsystem.DirTool('09_0*');
            for i = 1:length(dt.fqdns)
                this.testObj.SL_AFC(dt.fqdns{i}, dt.dns{i})
            end
            popd(pwd0)
        end
	end

 	methods (TestClassSetup)
		function setupAFC(this)
 			import mlafc.*;
            this.registry = AFCRegistry.instance('initialize');
            this.WORK = this.registry.workpath;
            setenv('WORK', this.WORK)
 			%this.testObj_ = KaysAFC;
 		end
	end

 	methods (TestMethodSetup)
		function setupAFCTest(this)
 			this.testObj = this.testObj_;
 			this.addTeardown(@this.cleanTestMethod);
 		end
    end
    
    %% PRIVATE

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanTestMethod(this)
        end
        function report_out(~, varargin)
            ip = inputParser;
            addRequired(ip, 'outmap', @(x) isa(x, 'containers.Map'))
            parse(ip, varargin{:})
            ipr = ip.Results;

            m = [];
            x = [];
            for pid = {'PT15' 'PT28' 'PT36'}
                m = [m; ipr.outmap(pid{1}).mean];                    
                x = [x; ipr.outmap(pid{1}).max];
            end
            fprintf('Engel Ia mean: %g\n', mean(m));
            fprintf('Engel Ia max:  %g\n', max(x));

            fprintf('Engel IIa mean: %g\n', ipr.outmap('PT34').mean);
            fprintf('Engel IIa max:  %g\n', ipr.outmap('PT34').max);                

            for pid = {'PT26' 'PT29' 'PT35'}
                m = [m; ipr.outmap(pid{1}).mean];                    
                x = [x; ipr.outmap(pid{1}).max];
            end
            fprintf('Engel IV mean: %g\n', mean(m));
            fprintf('Engel IV max:  %g\n', max(x));            
        end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

