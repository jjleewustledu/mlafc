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
        WORK = '/data/nil-bluearc/shimony/jjlee/FocalEpilepsy'
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
                [sims(ip),featifc,funcifc] = this.testObj.makeDice( ...
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
            setenv('WORK', this.WORK)
 			this.testObj_ = AFC;
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
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

