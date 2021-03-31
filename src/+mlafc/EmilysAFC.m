classdef EmilysAFC < mlafc.AFC
	%% EMILYSAFC  

	%  $Revision$
 	%  was created 10-Nov-2020 12:50:31 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlafc/src/+mlafc.
 	%% It was developed on Matlab 9.7.0.1471314 (R2019b) Update 7 for MACI64.  Copyright 2020 John Joowon Lee.
 	
	properties
        betFrac = 0.65 % for resection images only; smaller leaves more brain
        containerResections
        containerMprs % no resections
        containerT2w % no resections
        containerSegs
        containerAnatAves
        containerBolds
    end
    
    properties (Dependent)
        homeEmily
        homeGrid
        keys
    end
    
    methods (Static)
        function buildAfcProbOnly(varargin)
            ip = inputParser;
            addOptional(ip, 'toglob', 'PT*', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            atl = mlfourd.ImagingContext2( ...
                fullfile(getenv('REFDIR'), '711-2B_333.4dfp.hdr'));
            msk = mlfourd.ImagingContext2( ...
                '/Users/jjlee/Box/Leuthardt_Epilepsy_Project/Segmentations_06_and_07_2020/glm_atlas_mask_333.4dfp.hdr');
            for g = globFoldersT(ipr.toglob)
                pwd0 = pushd(fullfile(g{1}, 'AFC_and_resections', ''));

                afc = mlfourd.ImagingContext2('afc.4dfp.hdr');
                prob = copy(atl.fourdfp);
                prob.img = afc.fourdfp.img;
                prob.img(isnan(prob.img)) = 0;
                prob = mlfourd.ImagingContext2(prob);
                prob = (prob + 1)./2; % tanh to sigmoid
                prob = prob.ones - prob; % prob of aberrancy
                prob = prob .* msk;
                assert(dipmax(prob) <= 1)
                assert(dipmin(prob) >= 0)
                prob.filepath = pwd;
                prob.fileprefix = 'afc_prob';
                prob.save()
                prob.nifti.save()
                
                popd(pwd0)
            end
        end
        function buildAfcProbOnResection(varargin)
            ip = inputParser;
            addOptional(ip, 'toglob', 'PT*', @ischar)
            addParameter(ip, 'cost', 'corratio', @ischar)
            addParameter(ip, 'dof', 12, @isscalar)
            addParameter(ip, 'allbrain', false, @islogical)
            parse(ip, varargin{:})
            ipr = ip.Results;
            ipr.dof = num2str(ipr.dof);
            
            for g = globFoldersT(ipr.toglob)
                if strcmp(g{1}, 'PT26') || strcmp(g{1}, 'PT28')
                    warning('mlafc:RuntimeWarning', ...
                        'EmilysAFC.buildAfcProbOnResection:  manual registration needed for %s', g{1})
                    continue
                end
                
                pwd0 = pushd(g{1});
                
                opts = ['-bins 256 -cost ' ipr.cost ' -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof ' ipr.dof ' -interp trilinear'];
                opts1 = '-paddingsize 0.0 -interp trilinear';
                resection = globT([g{1} '_*_FLIRT_111.nii.gz']);
                assert(~isempty(resection))
                segmentation = globT([g{1} '_*_segmentation_final_111.nii.gz']);
                assert(~isempty(segmentation))
                
                try
                    % flirt BOLD on MPR
                    [~,r] = mlbash(sprintf('flirt -in %s%s -ref %s%s -out %s%s -omat %s%s %s', ...
                        g{1}, '_anat_ave_t88_333_brain.nii.gz', ...
                        g{1}, '_mpr1_on_TRIO_Y_NDC_111_brain.nii.gz', ...
                        g{1}, '_anat_ave_t88_333_brain_on_mpr1_on_111.nii.gz', ...
                        g{1}, '_anat_ave_t88_333_brain_on_mpr1_on_111.mat', ...
                        opts));
                    % flirt MPR on resection
                    if ipr.allbrain
                        [~,r] = mlbash(sprintf('flirt -in %s%s -ref %s -out %s%s -omat %s%s %s', ...
                            g{1}, '_mpr1_on_TRIO_Y_NDC_111_brain.nii.gz', ...
                            resection{1}, ...
                            g{1}, '_mpr1_on_FLIRT_111_brain.nii.gz', ...
                            g{1}, '_mpr1_on_FLIRT_111.mat', ...
                            opts));
                    else
                        [~,r] = mlbash(sprintf('flirt -in %s%s -ref %s -out %s%s -omat %s%s %s', ...
                            g{1}, '_mpr1_on_TRIO_Y_NDC_111.nii.gz', ...
                            resection{1}, ...
                            g{1}, '_mpr1_on_FLIRT_111.nii.gz', ...
                            g{1}, '_mpr1_on_FLIRT_111.mat', ...
                            opts));
                    end

                    % compose xfm
                    [~,r] = mlbash(sprintf('convert_xfm -omat %s%s -concat %s%s %s%s', ...
                        g{1}, '_anat_ave_on_FLIRT_111.mat', ...
                        g{1}, '_mpr1_on_FLIRT_111.mat ', ...
                        g{1}, '_anat_ave_t88_333_brain_on_mpr1_on_111.mat'));

                    % flirt afc_prob on resection
                    afc_prob_fn = 'AFC_and_resections/afc_prob.nii.gz';
                    assert(isfile(afc_prob_fn))
                    [~,r] = mlbash(sprintf('flirt -in %s -applyxfm -init %s%s -out %s%s %s -ref %s', ...
                        afc_prob_fn, ...
                        g{1}, '_anat_ave_on_FLIRT_111.mat', ...
                        g{1}, '_afc_prob_111.nii.gz ', ...
                        opts1, ...
                        resection{1}));
                    
                    qc = mlfourd.ImagingContext2([g{1} '_afc_prob_111.nii.gz']);
                    qc.fsleyes(resection{1}, segmentation{1})
                catch ME
                    warning('mlafc:RuntimeError', r)
                    handerror(ME)
                end
                
                popd(pwd0)
            end
        end
        function buildFCToResection(varargin)
            import mlafc.EmilysAFC.applyxfm
            
            ip = inputParser;
            addOptional(ip, 'toglob', 'PT*', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;            
            
            atl = mlfourd.ImagingContext2(fullfile(getenv('REFDIR'), '711-2B_333.nii.gz'));
            gm3d = mlfourd.ImagingContext2(fullfile(getenv('REFDIR'), 'gm3d_and_cerebellum.nii.gz'));
            gmbin = reshape(logical(gm3d), [48*64*48 1]);
            Ngm = dipsum(gmbin); % ~ 51856
            
            for g = globFoldersT(ipr.toglob)
                pwd0 = pushd(g{1});
                
                segfn = globT([g{1} '_*_segmentation_final_111.nii.gz']);
                seg = mlfourd.ImagingContext2(segfn{1});
                mat1to3 = fullfile(getenv('REFDIR'), '711-2B_111_on_333.mat');
                seg = applyxfm(seg, mat1to3, [g{1} '_segmentation_final_333'], atl);
                seg.fsleyes([g{1} '_anat_ave_t88_333.nii.gz'])
                segbin = logical(seg);
                segbin = flip(segbin, 2);
                segbin = reshape(segbin, [48*64*48 1]); 
                
                mat = load(fullfile('/Users/jjlee/Box/DeepNetFCProject/Epilepsy/NoiseInjection', ...
                                    g{1}, ...
                                    [g{1} '_BOLD.mat']));
                
                seed = mean(mat.dat(segbin,:), 1); % 1 x Nt
                
                dat = mat.dat(gmbin,:); % Ngm x Nt
                imgfc = zeros(Ngm, 1);
                imgpval = zeros(Ngm, 1);
                for p = 1:Ngm
                    [cc,pval] = corrcoef(seed', dat(p,:)'); % col vec, col vec
                    imgfc(p) = cc(1,2);
                    imgpval(p) = pval(1,2);
                end
                fc = copy(atl.zeros);
                fc.filepath = pwd;
                fc.fileprefix = [g{1} '_fc_to_resection'];
                fc = fc.nifti;
                fc.img(gmbin) = imgfc;
                fc.fsleyes()
                fc.save()
                pval = copy(atl.zeros);
                pval.filepath = pwd;
                pval.fileprefix = [g{1} '_pval_to_resection'];
                pval = pval.nifti;
                pval.img(gmbin) = imgpval;
                pval.fsleyes()
                pval.save()
                
                popd(pwd0)
            end
        end
        function buildSoftmaxOnResection(varargin)
            ip = inputParser;
            addOptional(ip, 'toglob', 'PT*', @ischar)
            addParameter(ip, 'cost', 'corratio', @ischar)
            addParameter(ip, 'dof', 12, @isscalar)
            addParameter(ip, 'search', 45, @isscalar)
            addParameter(ip, 'finesearch', 10, @isscalar)
            addParameter(ip, 'allbrain', true, @islogical)
            parse(ip, varargin{:})
            ipr = ip.Results;
            ipr.dof = num2str(ipr.dof);
            
            for g = globFoldersT(ipr.toglob)  
                if strcmp(g{1}, 'PT26') || strcmp(g{1}, 'PT28')
                    warning('mlafc:RuntimeWarning', ...
                        'EmilysAFC.buildSoftmaxOnResection:  manual registration needed for %s', g{1})
                    continue
                end
                
                pwd0 = pushd(g{1});
                
                % ensure nii.gz
                if ~isfile('AFC_product.nii.gz')
                    assert(isfile('AFC_product.4dfp.hdr'))
                    ic = mlfourd.ImagingContext2('AFC_product.4dfp.hdr');
                    ic.nifti.save()
                end
                
                s = sprintf('-%i %i', ipr.search, ipr.search);
                fs = sprintf('%i', ipr.finesearch);
                opts = ['-bins 256 -cost ' ipr.cost ' -searchrx ' s ' -searchry ' s ' -searchrz ' s ' -dof ' ipr.dof ' -finesearch ' fs ' -interp trilinear'];
                opts_6 = ['-bins 256 -cost ' ipr.cost ' -searchrx ' s ' -searchry ' s ' -searchrz ' s ' -dof 6 -interp trilinear']; 
                opts_6_fine = ['-bins 256 -cost ' ipr.cost ' -searchrx ' s ' -searchry ' s ' -searchrz ' s ' -dof 6 -finesearch 5 -interp trilinear']; 
                opts_12 = ['-bins 256 -cost ' ipr.cost ' -searchrx ' s ' -searchry ' s ' -searchrz ' s ' -dof 12 -interp trilinear'];
                opts_trans = ['-bins 256 -cost ' ipr.cost ' -searchrx ' s ' -searchry ' s ' -searchrz ' s ' -dof 6 -schedule /usr/local/fsl/etc/flirtsch/sch3Dtrans_3dof -interp trilinear'];
                opts_apply = '-paddingsize 0.0 -interp trilinear';
                resection = globT([g{1} '_*_FLIRT_111.nii.gz']);
                assert(~isempty(resection))
                segmentation = globT([g{1} '_*_segmentation_final_111.nii.gz']);
                assert(~isempty(segmentation))
                
                try
                    % flirt BOLD on MPR
                    [~,r] = mlbash(sprintf('flirt -in %s%s -ref %s%s -out %s%s -omat %s%s %s', ...
                        g{1}, '_anat_ave_t88_333_brain.nii.gz', ...
                        g{1}, '_mpr1_on_TRIO_Y_NDC_111_brain.nii.gz', ...
                        g{1}, '_anat_ave_t88_333_brain_on_mpr1_on_111.nii.gz', ...
                        g{1}, '_anat_ave_t88_333_brain_on_mpr1_on_111.mat', ...
                        opts_12));
                    % flirt MPR on resection
                    if ipr.allbrain
                        [~,r] = mlbash(sprintf('flirt -in %s%s -ref %s -out %s%s -omat %s%s %s', ...
                            g{1}, '_mpr1_on_TRIO_Y_NDC_111_brain.nii.gz', ...
                            resection{1}, ...
                            g{1}, '_mpr1_on_FLIRT_111_brain.nii.gz', ...
                            g{1}, '_mpr1_on_FLIRT_111.mat', ...
                            opts));
                    else
                        [~,r] = mlbash(sprintf('flirt -in %s%s -ref %s -out %s%s -omat %s%s %s', ...
                            g{1}, '_mpr1_on_TRIO_Y_NDC_111.nii.gz', ...
                            resection{1}, ...
                            g{1}, '_mpr1_on_FLIRT_111.nii.gz', ...
                            g{1}, '_mpr1_on_FLIRT_111.mat', ...
                            opts));
                    end

                    % compose xfm
                    [~,r] = mlbash(sprintf('convert_xfm -omat %s%s -concat %s%s %s%s', ...
                        g{1}, '_anat_ave_on_FLIRT_111.mat', ...
                        g{1}, '_mpr1_on_FLIRT_111.mat ', ...
                        g{1}, '_anat_ave_t88_333_brain_on_mpr1_on_111.mat'));

                    % flirt afc_prob on resection
                    afc_prob_fn = 'AFC_product.nii.gz';
                    assert(isfile(afc_prob_fn))
                    [~,r] = mlbash(sprintf('flirt -in %s -applyxfm -init %s%s -out %s%s %s -ref %s', ...
                        afc_prob_fn, ...
                        g{1}, '_anat_ave_on_FLIRT_111.mat', ...
                        g{1}, '_afc_prob_111.nii.gz ', ...
                        opts_apply, ...
                        resection{1}));
                    
                    qc = mlfourd.ImagingContext2([g{1} '_afc_prob_111.nii.gz']);
                    qc.fsleyes(resection{1}, segmentation{1})
                catch ME
                    warning('mlafc:RuntimeError', r)
                    handerror(ME)
                end
                
                popd(pwd0)
            end
        end         
        function calc_dice(varargin)
            ip = inputParser;
            addOptional(ip, 'toglob', 'PT*', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            for g = globFoldersT(ipr.toglob)
                pwd0 = pushd(g{1});
                
                segmentation = globT([g{1} '_*_segmentation_final_111.nii.gz']);
                seg = mlfourd.ImagingContext2(segmentation{1});
                afc_prob = mlfourd.ImagingContext2([g{1} '_afc_prob_111.nii.gz']);
                gm = mlfourd.ImagingContext2(fullfile(getenv('REFDIR'), 'gm3d_111.nii.gz')); % no cerebellum
                gm = gm.masked(double(afc_prob.numgt(0.008))); % 0.008 is the left tail of histograms
%                thr = 0.00891826; % mean of modes for Engel I, II
%                afc_prob = afc_prob.thresh(thr);
                d = afc_prob.dice(seg, gm);
                d = d.nifti.img;
                fprintf('dice(%s) = %g\n', afc_prob.fileprefix, d)
                pr = afc_prob.nifti.img(logical(gm));
                fprintf('mode(%s) = %g\n', afc_prob.fileprefix, mode(pr))
                
                popd(pwd0)
            end
        end        
        function calc_histogram(varargin)
            ip = inputParser;
            addOptional(ip, 'toglob', 'PT*', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            for g = globFoldersT(ipr.toglob)
                pwd0 = pushd(g{1});
                
                afc_prob = mlfourd.ImagingContext2([g{1} '_afc_prob_111.nii.gz']);                
                afc_prob = afc_prob.thresh(0.001);
                h = figure;
                histo = histogram(afc_prob, 'DisplayStyle', 'stairs', 'LineWidth', 2);
                ax = gca;
                ax.XRuler.Exponent = 0;
                xtickformat('%.3f')
                set(gca, 'fontsize', 12)
                xlabel('probability dissimilarity', 'FontSize', 18)
                ylabel('number voxels', 'FontSize', 18)
                m = mode(afc_prob.nifti.img(logical(afc_prob)));
                s = sprintf('mode = %.4f\nfull width half max = %.4f', m, fwhm(histo)); 
                annotation('textbox', [0.16 .2 .5 0.1], 'String', s, 'FitBoxToText', 'on', 'FontSize', 16, 'LineStyle', 'none')
                
                savefig(h, ...
                    sprintf('EmilysAFC_calc_histogram_%s_afc_prob_111.fig', g{1}))
                figs = get(0, 'children');
                saveas(figs(1), ...
                    sprintf('EmilysAFC_calc_histogram_%s_afc_prob_111.png', g{1}))
                close(figs(1))
                
                popd(pwd0)
            end

            function w = fwhm(histo)
                hh = 0.5 * (max(histo.Values) - min(histo.Values));
                idx1 = find(histo.Values >= hh, 1, 'first');
                idx2 = find(histo.Values >= hh, 1, 'last');
                w = histo.BinEdges(idx2) - histo.BinEdges(idx1);
            end
        end
        function calc_jsdiv(varargin)
            ip = inputParser;
            addOptional(ip, 'toglob', 'PT*', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            for g = globFoldersT(ipr.toglob)
                pwd0 = pushd(g{1});
                
                segmentation = globT([g{1} '_*_segmentation_final_111.nii.gz']);
                seg = mlfourd.ImagingContext2(segmentation{1});
                afc_prob = mlfourd.ImagingContext2([g{1} '_afc_prob_111.nii.gz']);
                gm = mlfourd.ImagingContext2(fullfile(getenv('REFDIR'), 'gm3d_111.nii.gz')); % no cerebellum
%                gm = gm.masked(double(afc_prob.numgt(0.008))); % 0.008 is the left tail of histograms                
                jsdiv = afc_prob.jsdiv(seg, gm);
                fprintf('jsdiv(%s) = %g\n',afc_prob.fileprefix, jsdiv)
                
                popd(pwd0)
            end
        end
        function calc_perfcurve(varargin)
            ip = inputParser;
            addOptional(ip, 'toglob', 'PT*', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            for g = globFoldersT(ipr.toglob)
                pwd0 = pushd(g{1});
                
                segmentation = globT([g{1} '_*_segmentation_final_111.nii.gz']);
                seg = mlfourd.ImagingContext2(segmentation{1});
                afc_prob = mlfourd.ImagingContext2([g{1} '_afc_prob_111.nii.gz']);
                gm = mlfourd.ImagingContext2(fullfile(getenv('REFDIR'), 'gm3d_111.nii.gz')); % no cerebellum
                msk = logical(gm);
                                
                seg_ = double(seg.nifti.img(msk));
                afc_prob_ = double(afc_prob.nifti.img(msk));
                [x,y,t,auc,OPTROCPT] = perfcurve(seg_, afc_prob_, 1);
                plot(x,y)
                hold on
                plot(OPTROCPT(1),OPTROCPT(2),'ro')
                xlabel('False positive rate') 
                ylabel('True positive rate')
                title(sprintf('ROC Curve for %s\nAUC = %g', ...
                              afc_prob.fileprefix, auc))
                hold off
                
                popd(pwd0)
            end
        end
        function calc_var(varargin)
            ip = inputParser;
            addOptional(ip, 'toglob', 'PT*', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            for g = globFoldersT(ipr.toglob)
                pwd0 = pushd(g{1});
                
                afc_prob = mlfourd.ImagingContext2([g{1} '_afc_prob_111.nii.gz']);
                v = var(afc_prob.nifti.img(logical(afc_prob.numgt(0.001))));
                fprintf('var(%s) = %g\n',afc_prob.fileprefix, v)
                fprintf('std(%s) = %g\n',afc_prob.fileprefix, sqrt(v))
                
                popd(pwd0)
            end
        end
        function [sim,featic,funcic] = calcdice__(featdata, funcdata, varargin)
            
            import mlafc.EmilysAFC.flip12
            
            ip = inputParser;
            addRequired(ip, 'featdata', @isnumeric)
            addRequired(ip, 'funcdata', @isnumeric)   
            addParameter(ip, 'Nsigma', 2, @isnumeric)
            parse(ip, featdata, funcdata, varargin{:})
            ipr = ip.Results;
            
            [GLMmask,~,~,glmic] = mlperceptron.PerceptronRegistry.read_glm_atlas_mask();
            GLMmask(find(GLMmask)) = 1; 
            
            % reshape
            GLMmask  = reshape(GLMmask,  [48 64 48]);
            featdata = reshape(featdata, [48 64 48]) .* GLMmask; 
            funcdata = reshape(funcdata, [48 64 48]) .* GLMmask; 
            
            % threshold funcdata
            mu = nanmean(funcdata(logical(GLMmask)));
            sigma = nanstd(funcdata(logical(GLMmask)));
            
            funcdata = funcdata < mu - ipr.Nsigma*sigma;
            
            featifc = copy(glmic.fourdfp);
            featifc.img = flip12(single(featdata .* GLMmask));
            featifc.filepath = pwd;
            featifc.fileprefix = 'calcdice_label';
            featic = mlfourd.ImagingContext2(featifc);
            funcifc = copy(glmic.fourdfp);
            funcifc.img = flip12(single(funcdata .* GLMmask));
            funcifc.filepath = pwd;
            funcifc.fileprefix = sprintf('calcdice_afc_%isd', ipr.Nsigma);
            funcic = mlfourd.ImagingContext2(funcifc);
            
            sim = dice(logical(featdata), logical(funcdata));  
            fprintf('mlafc.AFC.calcdice.sim -> %g\n', sim)
        end
        function img = flip12(img)
            img = flip(flip(img, 1), 2);
        end
        function visualizeAfcProb(varargin)
            ip = inputParser;
            addOptional(ip, 'toglob', 'PT*', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            for g = globFoldersT(ipr.toglob)
                pwd0 = pushd(g{1});
                
                resection = globT([g{1} '_*_FLIRT_111.nii.gz']);
                assert(~isempty(resection))
                segmentation = globT([g{1} '_*_segmentation_final_111.nii.gz']);
                assert(~isempty(segmentation))
                edge_seg = mlfourd.ImagingFormatContext(segmentation{1});
                edge_seg.img = edge3(edge_seg.img, 'approxcanny', 0.2);
                edge_seg.fileprefix = [strrep(segmentation{1}, '.nii.gz', '') '_edge'];
%                edge_seg.save()
                
                qc = mlfourd.ImagingContext2([g{1} '_afc_prob_111.nii.gz']);
                qc.fsleyes(resection{1}, edge_seg.filename)
                
                popd(pwd0)
            end
        end
        
        %% UTILITIES
        
        function ic = bet(ic0, betFrac)
            bin = fullfile(getenv('FSLDIR'), 'bin', 'bet');
            if contains(ic0.fileprefix, 'mpr')
                t2_fqfp = strrep(ic0.fqfp, 'mpr1', 't2w');
                assert(isfile([t2_fqfp '.nii.gz']))
                mlbash(sprintf('%s %s %s_brain -A2 %s -f %g -g 0 -m', bin, ic0.fqfp, ic0.fqfp, t2_fqfp, betFrac))
                BET = fullfile(ic0.filepath, 'BET', '');
                ensuredir(BET)
                movefile(fullfile(ic0.filepath, '*_mask.*'), BET, 'f')
                movefile(fullfile(ic0.filepath, '*_mesh.*'), BET, 'f')
            else
                mlbash(sprintf('%s %s %s_brain -R -f %i -g 0 -m', bin, ic0.fqfp, ic0.fqfp, betFrac))
            end
            ic = mlfourd.ImagingContext2([ic0.fqfp '_brain.nii.gz']);
            ic.selectImagingFormatTool()
        end
        function ic = betZ(ic0, betFrac)
            bin = fullfile(getenv('FSLDIR'), 'bin', 'bet');
            mlbash(sprintf('%s %s %s_brain -Z -f %i -g 0 -m', bin, ic0.fqfp, ic0.fqfp, betFrac))
            ic = mlfourd.ImagingContext2([ic0.fqfp '_brain.nii.gz']);
            ic.selectImagingFormatTool()
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
        function ic = applyxfm(ic0, mat, out_fqfp, icref)
            bin = fullfile(getenv('FSLDIR'), 'bin', 'flirt');
            mlbash(sprintf('%s -in %s -applyxfm -init %s -out %s -paddingsize 0.0 -interp nearestneighbour -ref %s', bin, ic0.fqfp, mat, out_fqfp, icref.fqfp))
            %flirt -in PT15_4_seg_111_nopriorsurg.nii.gz -applyxfm -init PT15_4_FLIRT_on_mpr1_111_brain.mat -out PT15_4_seg_111_nopriorsurg_on_mpr_111.nii.gz -paddingsize 0.0 -interp nearestneighbour -ref PT15_mpr1_on_TRIO_Y_NDC_111.nii.gz
            ic = mlfourd.ImagingContext2([out_fqfp '.nii.gz']);
            ic.selectImagingFormatTool()
        end
        function [icres,icseg] = antsToMprDeformed(icres0, icseg0)
            %% see also:
            %  https://github.com/ANTsX/ANTs/wiki/Forward-and-inverse-warps-for-warping-images,-pointsets-and-Jacobians
            
            [~,key] = fileparts(pwd);
            syn = fullfile(getenv('ANTSPATH'), 'antsRegistrationSyNQuick.sh');
            apply = fullfile(getenv('ANTSPATH'), 'antsApplyTransforms');
            resfn = [key '_resection_toMprDeformed.nii.gz'];
            segfn = [key '_seg_toMprDeformed.nii.gz'];
            icmpr = mlfourd.ImagingContext2([key '_mpr1_on_TRIO_Y_NDC_111_brain.nii.gz']);
            assert(isfile(icmpr.fqfn))
            
            mlbash(sprintf('%s -d 3 -f %s -m %s -o movingToMpr_ -t s', syn, icmpr.fqfn, icres0.fqfn))
            mlbash(sprintf('%s -d 3 -i %s -r %s -t movingToMpr_1Warp.nii.gz -t movingToMpr_0GenericAffine.mat -o %s', apply, icres0.fqfn, icmpr.fqfn, resfn))
            mlbash(sprintf('%s -d 3 -i %s -r %s -t movingToMpr_1Warp.nii.gz -t movingToMpr_0GenericAffine.mat -n GenericLabel -o %s', apply, icseg0.fqfn, icmpr.fqfn, segfn))
            %antsRegistrationSyNQuick.sh -d 3 -f PT15_mpr1_on_TRIO_Y_NDC_111.nii.gz -m PT15_anat_ave_t88_333.nii.gz -o movingToFixed_ -t s
            %antsApplyTransforms -d 3 -i PT15_resection_on_mpr_111.nii.gz -r PT15_anat_ave_t88_333.nii.gz -t [movingToFixed_0GenericAffine.mat, 1] -t movingToFixed_1InverseWarp.nii.gz -o PT15_4_FLIRT_toAnatDeformed.nii.gz
            %antsApplyTransforms -d 3 -i PT15_seg_on_mpr_111.nii.gz -r PT15_anat_ave_t88_333.nii.gz -t [movingToFixed_0GenericAffine.mat, 1] -t movingToFixed_1InverseWarp.nii.gz -n GenericLabel -o PT15_4_seg_toAnatDeformed.nii.gz
            icres = mlfourd.ImagingContext2(resfn);
            icres.selectImagingFormatTool()
            icseg = mlfourd.ImagingContext2(segfn);
            icseg.selectImagingFormatTool()
        end
        function [icres,icseg] = antsToBoldDeformed(icres0, icseg0, icmpr, icbold)
            %% see also:
            %  https://github.com/ANTsX/ANTs/wiki/Forward-and-inverse-warps-for-warping-images,-pointsets-and-Jacobians
            
            [~,key] = fileparts(pwd);
            syn = fullfile(getenv('ANTSPATH'), 'antsRegistrationSyNQuick.sh');
            apply = fullfile(getenv('ANTSPATH'), 'antsApplyTransforms');
            resfn = [key '_resection_toBoldDeformed.nii.gz'];
            segfn = [key '_seg_toBoldDeformed.nii.gz'];
            
            mlbash(sprintf('%s -d 3 -f %s -m %s -o movingToMpr2_ -t s', syn, icmpr.fqfn, icbold.fqfn))
            mlbash(sprintf('%s -d 3 -i %s -r %s -t [movingToMpr2_0GenericAffine.mat, 1] -t movingToMpr2_1InverseWarp.nii.gz -o %s', apply, icres0.fqfn, icbold.fqfn, resfn))
            mlbash(sprintf('%s -d 3 -i %s -r %s -t [movingToMpr2_0GenericAffine.mat, 1] -t movingToMpr2_1InverseWarp.nii.gz -n GenericLabel -o %s', apply, icseg0.fqfn, icbold.fqfn, segfn))
            %antsRegistrationSyNQuick.sh -d 3 -f PT15_mpr1_on_TRIO_Y_NDC_111.nii.gz -m PT15_anat_ave_t88_333.nii.gz -o movingToFixed_ -t s
            %antsApplyTransforms -d 3 -i PT15_resection_on_mpr_111.nii.gz -r PT15_anat_ave_t88_333.nii.gz -t [movingToFixed_0GenericAffine.mat, 1] -t movingToFixed_1InverseWarp.nii.gz -o PT15_4_FLIRT_toAnatDeformed.nii.gz
            %antsApplyTransforms -d 3 -i PT15_seg_on_mpr_111.nii.gz -r PT15_anat_ave_t88_333.nii.gz -t [movingToFixed_0GenericAffine.mat, 1] -t movingToFixed_1InverseWarp.nii.gz -n GenericLabel -o PT15_4_seg_toAnatDeformed.nii.gz
            icres = mlfourd.ImagingContext2(resfn);
            icres.selectImagingFormatTool();
            icseg = mlfourd.ImagingContext2(segfn);
            icseg.selectImagingFormatTool();
        end
    end

	methods 
        
        %% GET
        
        function g = get.homeEmily(this)
            if contains(hostname, 'precuneal')
                g = '/Users/jjlee/Box/Leuthardt_Epilepsy_Project';
                return
            end
            g = '/data/nil-bluearc/shimony/jjlee/FocalEpilepsy/Emily';
        end
        function g = get.homeGrid(this)
            if contains(hostname, 'precuneal')
                g = '/Users/jjlee/Box/DeepNetFCProject/Epilepsy/NoiseInjection';
                return
            end
            g = '/data/nil-bluearc/shimony/jjlee/DeepNetFCProject/Epilepsy';         
        end
        function g = get.keys(this)
            g = this.containerResections.keys;
        end
        
        %% 
        
 		function this = EmilysAFC(varargin)
 			this = this@mlafc.AFC(varargin{:});
            
            this.containerResections = this.buildContainerResections();
            this.containerMprs = this.buildContainerMprs();
            this.containerT2w = this.buildContainerT2w();
            this.containerSegs = this.buildContainerSegs();
            this.containerBolds = this.buildContainerBolds();
            this.containerAnatAves = this.buildContainerAnatAves();
        end
        
        function [res,seg] = buildSegOnBold(this, key)
            %% uses bet, flirt, antsRegistrationSyNQuick.sh to register resection & segmentation on BOLD 
            
            import mlfourd.ImagingContext2
            assert(isfolder(key))
            
            pwd0 = pushd(key);
            if contains(key, filesep)
                [~,key] = fileparts(key);
            end
            res = ImagingContext2(this.containerResections(key));
            mpr = ImagingContext2(this.containerMprs(key));
            seg = ImagingContext2(this.containerSegs(key));
            bold = ImagingContext2([key '_anat_ave_t88_333.nii.gz']);
            res = this.bet(res, this.betFrac);
            mpr = this.bet(mpr, 0.5);
            bold = this.bet(bold, 0.5);
            [res,mat9] = this.flirt(9, res, mpr, sprintf('%s_resection_on_mpr_111', key));
            seg = this.applyxfm(seg, mat9, sprintf('%s_seg_on_mpr_111', key), mpr);
            %[res,mat6] = this.flirt(6, res, bold, sprintf('%s_resection_on_anat_333', key));
            %seg = this.applyxfm(seg, mat6, sprintf('%s_seg_on_anat_333', key), bold);
            %[res,seg] = this.antsToMprDeformed(res, seg);
            [res,seg] = this.antsToBoldDeformed(res, seg, mpr, bold);
            popd(pwd0)
        end
        function [res,seg] = buildSegOnBold1(this, key)
            %% uses bet, flirt, antsRegistrationSyNQuick.sh to register resection & segmentation on BOLD 
            
            import mlfourd.ImagingContext2
            assert(isfolder(key))
            
            pwd0 = pushd(key);
            if contains(key, filesep)
                [~,key] = fileparts(key);
            end
            res = ImagingContext2(this.containerResections(key));
            mpr = ImagingContext2(this.containerMprs(key));
            seg = ImagingContext2(this.containerSegs(key));
            bold = ImagingContext2([key '_anat_ave_t88_333.nii.gz']);
            res = this.betZ(res, this.betFrac);
            mpr = this.bet(mpr, 0.5);
            bold = this.bet(bold, 0.5);
            [res,mat9] = this.flirt(9, res, mpr, sprintf('%s_resection_on_mpr_111', key));
            seg = this.applyxfm(seg, mat9, sprintf('%s_seg_on_mpr_111', key), mpr);
            %[res,mat6] = this.flirt(6, res, bold, sprintf('%s_resection_on_anat_333', key));
            %seg = this.applyxfm(seg, mat6, sprintf('%s_seg_on_anat_333', key), bold);
            %[res,seg] = this.antsToMprDeformed(res, seg);
            [res,seg] = this.antsToBoldDeformed(res, seg, mpr, bold);
            popd(pwd0)
        end
        function [res,seg] = buildSegOnBold_rigid(this, key)
            %% uses bet, flirt, antsRegistrationSyNQuick.sh to register resection & segmentation on BOLD 
            
            import mlfourd.ImagingContext2
            assert(isfolder(key))
            
            pwd0 = pushd(key);
            if contains(key, filesep)
                [~,key] = fileparts(key);
            end
            res = ImagingContext2(this.containerResections(key));
            mpr = ImagingContext2(this.containerMprs(key));
            seg = ImagingContext2(this.containerSegs(key));
            bold = ImagingContext2([key '_anat_ave_t88_333.nii.gz']);
            %res = this.bet(res, this.betFrac);
            %mpr = this.bet(mpr, 0.5);
            %bold = this.bet(bold, 0.7);
            [res,mat6] = this.flirt(6, res, mpr, sprintf('%s_resection_on_mpr_111', key));
            seg = this.applyxfm(seg, mat6, sprintf('%s_seg_on_mpr_111', key), mpr);
            %[res,mat6] = this.flirt(6, res, bold, sprintf('%s_resection_on_anat_333', key));
            %seg = this.applyxfm(seg, mat6, sprintf('%s_seg_on_anat_333', key), bold);
            %[res,seg] = this.antsToMprDeformed(res, seg);
            [res,seg] = this.antsToBoldDeformed(res, seg, mpr, bold);
            popd(pwd0)
        end
        function [sim,featic,funcic] = makeDice(this, varargin)
            %% MAKEDICE
            %  Preconditions:
            %  -- completion of SL_fMRI_initialization
            %  -- completion of mlperceptron.PerceptronRelease factory
            %  @param patientdir is folder.
            %  @param patientid is char.
            %  @param afc_filename is file.
            %  @param feature_filename is file.
            %  @returns similarity
            %  @returns ImagingFormatContext
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired( ip, 'patientdir', @isfolder)
            addRequired( ip, 'patientid', @ischar)
            addParameter(ip, 'afc_filename', ['AFC_this_' datestr(now, 30) '.mat'], @isfile)
            addParameter(ip, 'feature_filename', '', @isfile)
            addParameter(ip, 'load_afc', true, @islogical)
            addParameter(ip, 'Nsigma', 2, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.patientdir = ipr.patientdir;
            this.patientid = ipr.patientid;
            
            if ipr.load_afc && isfile(ipr.afc_filename)
                mat = load(ipr.afc_filename');
                this.afc_map = mat.obj.afc_map;
            else
                warning('mlafc:RuntimeWarning', 'makeDice could not file afc file, so computing with makeSearchlightMap')
                this = this.makeSearchlightMap(varargin{:});
            end
            
            this.feature = this.imgread(ipr.feature_filename);
            [sim,featic,funcic] = this.calcdice(single(this.feature), ...
                                                single(this.afc_map), ...
                                                'Nsigma', ipr.Nsigma);
            featic.fourdfp.save()
            funcic.fourdfp.save()
            featic.nifti.save()
            funcic.nifti.save()
        end
    end 
    
    %% PRIVATE
    
    methods (Access = private, Static)
        function ic = rescueOriginator111(varargin)
            ic = mlfourd.ImagingContext2(varargin{:});
            nii = ic.nifti;
            nii.originator = [88 104 88];
            nii.save()
            ic = mlfourd.ImagingContext2(nii);
        end
        function ic = rescueOriginator333(varargin)
            ic = mlfourd.ImagingContext2(varargin{:});
            nii = ic.nifti;
            nii.originator = [72 96 72];
            nii.save()
            ic = mlfourd.ImagingContext2(nii);
        end
        function obj = rescue_AFC_and_resections_mat(fn)
            mat = load(fn);
            obj.ipr = mat.ipr;
            for p = properties(mat.this)'
                obj.(p{1}) = mat.this.(p{1});
            end
            [pth,fp] = fileparts(fn);
            fn1 = fullfile(pth, [fp '_obj.mat']);
            save(fn1, 'obj')
        end        
    end
    
    methods (Access = private)
        function c = buildContainerResections(this)
            %% build Map of FLAIRs of resections from final surgery.
            %  @return Map of ImagingContext2.
            
            c = containers.Map;
            import mlfourd.ImagingContext2
            home = fullfile(this.homeEmily, 'Segmentations_06_and_07_2020');
            c('PT15') = ImagingContext2(fullfile(home, 'PT15', 'PT15_4_FLIRT_111.nii.gz'));
            c('PT26') = ImagingContext2(fullfile(home, 'PT26', 'PT26_3_FLIRT_111.nii.gz'));
            c('PT28') = ImagingContext2(fullfile(home, 'PT28', 'PT28_4_FLIRT_111.nii.gz'));
            c('PT34') = ImagingContext2(fullfile(home, 'PT34', 'PT34_15_FLIRT_111.nii.gz'));
            c('PT35') = ImagingContext2(fullfile(home, 'PT35', 'PT35_6_FLIRT_111.nii.gz'));
            c('PT36') = ImagingContext2(fullfile(home, 'PT36', 'PT36_5_FLIRT_111.nii.gz')); 
        end
        function c = buildContainerMprs(this)
            %% build Map of MPRAGEs before final surgery.
            %  @return Map of ImagingContext2.
            
            c = containers.Map;
            import mlfourd.ImagingContext2
            home = fullfile(this.homeEmily, 'Segmentations_06_and_07_2020');
            c('PT15') = ImagingContext2(fullfile(home, 'PT15', 'PT15_mpr1_on_TRIO_Y_NDC_111.nii.gz'));
            c('PT26') = ImagingContext2(fullfile(home, 'PT26', 'PT26_mpr1_on_TRIO_Y_NDC_111.nii.gz'));
            c('PT28') = ImagingContext2(fullfile(home, 'PT28', 'PT28_mpr1_on_TRIO_Y_NDC_111.nii.gz'));
            c('PT34') = ImagingContext2(fullfile(home, 'PT34', 'PT34_mpr1_on_TRIO_Y_NDC_111.nii.gz'));
            c('PT35') = ImagingContext2(fullfile(home, 'PT35', 'PT35_mpr1_on_TRIO_Y_NDC_111.nii.gz'));
            c('PT36') = ImagingContext2(fullfile(home, 'PT36', 'PT36_mpr1_on_TRIO_Y_NDC_111.nii.gz')); 
        end
        function c = buildContainerT2w(this)
            %% build Map of T2w before final surgery.
            %  @return Map of ImagingContext2.
            
            c = containers.Map;
            import mlfourd.ImagingContext2
            home = fullfile(this.homeEmily, 'Segmentations_06_and_07_2020');
            c('PT15') = ImagingContext2(fullfile(home, 'PT15', 'PT15_t2w_on_TRIO_Y_NDC_111.nii.gz'));
            c('PT26') = ImagingContext2(fullfile(home, 'PT26', 'PT26_t2w_on_TRIO_Y_NDC_111.nii.gz'));
            c('PT28') = ImagingContext2(fullfile(home, 'PT28', 'PT28_t2w_on_TRIO_Y_NDC_111.nii.gz'));
            c('PT34') = ImagingContext2(fullfile(home, 'PT34', 'PT34_t2w_on_TRIO_Y_NDC_111.nii.gz'));
            c('PT35') = ImagingContext2(fullfile(home, 'PT35', 'PT35_t2w_on_TRIO_Y_NDC_111.nii.gz'));
            c('PT36') = ImagingContext2(fullfile(home, 'PT36', 'PT36_t2w_on_TRIO_Y_NDC_111.nii.gz')); 
        end
        function c = buildContainerSegs(this)
            %% build Map of labelled segmentations of resections from final surgery, excluding prior surgical resections.
            %  @return Map of ImagingContext2.
            
            c = containers.Map;
            import mlfourd.ImagingContext2
            home = fullfile(this.homeEmily, 'Segmentations_06_and_07_2020');
            c('PT15') = ImagingContext2(fullfile(home, 'PT15', 'PT15_4_seg_111_nopriorsurg.nii.gz'));
            c('PT26') = ImagingContext2(fullfile(home, 'PT26', 'PT26_3_segmentation_111.nii.gz'));
            c('PT28') = ImagingContext2(fullfile(home, 'PT28', 'PT28_4_seg_111_nopriorsurg.nii.gz'));
            c('PT34') = ImagingContext2(fullfile(home, 'PT34', 'PT34_15_Segmentation_111.nii.gz'));
            c('PT35') = ImagingContext2(fullfile(home, 'PT35', 'PT35_6_seg_111.nii.gz'));
            c('PT36') = ImagingContext2(fullfile(home, 'PT36', 'PT36_5_segmentation_111.nii.gz')); 
        end
        function c = buildContainerBolds(this)
            %% build Map of all of patient's BOLD.
            %  @return Map of mat filenames.
            
            c = containers.Map;
            home = this.homeGrid;
            c('PT15') = fullfile(home, 'PT15', 'PT15_BOLD.mat');
            c('PT26') = fullfile(home, 'PT26', 'PT26_BOLD.mat');
            c('PT28') = fullfile(home, 'PT28', 'PT28_BOLD.mat');
            c('PT34') = fullfile(home, 'PT34', 'PT34_BOLD.mat');
            c('PT35') = fullfile(home, 'PT35', 'PT35_BOLD.mat');
            c('PT36') = fullfile(home, 'PT36', 'PT36_BOLD.mat'); 
        end
        function c = buildContainerAnatAves(this)
            %% build Map of _anat_ave before final surgery.
            %  @return Map of ImagingContext2.
            
            c = containers.Map;
            import mlfourd.ImagingContext2
            home = fullfile(this.homeEmily, 'Segmentations_06_and_07_2020');
            c('PT15') = ImagingContext2(fullfile(home, 'PT15', 'PT15_anat_ave_t88_333.nii.gz'));
            c('PT26') = ImagingContext2(fullfile(home, 'PT26', 'PT26_anat_ave_t88_333.nii.gz'));
            c('PT28') = ImagingContext2(fullfile(home, 'PT28', 'PT28_anat_ave_t88_333.nii.gz'));
            c('PT34') = ImagingContext2(fullfile(home, 'PT34', 'PT34_anat_ave_t88_333.nii.gz'));
            c('PT35') = ImagingContext2(fullfile(home, 'PT35', 'PT35_anat_ave_t88_333.nii.gz'));
            c('PT36') = ImagingContext2(fullfile(home, 'PT36', 'PT36_anat_ave_t88_333.nii.gz')); 
        end        
        function ics = rescueOriginators(this, key)
            import mlfourd.ImagingContext2
            assert(isfolder(key))
            
            pwd0 = pushd(key);
            for item = globT('*_111*.nii.gz')
                this.rescueOriginator111(item{1})
            end
            for item = globT('*_333*.nii.gz')
                this.rescueOriginator333(item{1})
            end
            popd(pwd0)
        end
        
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

