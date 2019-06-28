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
        WORK = '/data/nil-bluearc/shimony/jjlee/Docker/SearchLight'
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
            this.testObj.SL_AFC( ...
                fullfile(getenv('WORK'), 'PT15', ''), 'PT15')
        end
        function test_New_Epilepsy(this)
            newepil = '/data/nil-bluearc/shimony/jjlee/New_Epilepsy';
            pwd0 = pushd(newepil);
            dt = mlsystem.DirTool('PT*');
            for i = 1:length(dt.fqdns)
                this.testObj.SL_AFC(dt.fqdns{i}, dt.dns{i})
            end
            this.testObj.SL_AFC(fullfile(newepil, '09_0009'), '09_0009')
            popd(pwd0)
        end
        function test_09_0009(this)
            newepil = '/data/nil-bluearc/shimony/jjlee/New_Epilepsy';
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

