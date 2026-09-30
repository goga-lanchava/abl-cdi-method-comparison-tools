classdef WalkthroughTest < matlab.unittest.TestCase
    % Regression tests for the figures quoted in WALKTHROUGH.md.
    %
    % The walkthrough commits to exact statistics for two example datasets,
    % which makes it a specification. These tests assert those values against
    % the example data in examples/, so the software and the article cannot
    % drift apart unnoticed.
    %
    % Run with:   runtests('tests')

    properties
        RepoRoot
        Examples
        AblFile
        CdiFile2026027
        CdiFile2026007
    end

    methods (TestClassSetup)
        function setupPaths(tc)
            here = fileparts(which('WalkthroughTest'));
            tc.RepoRoot = fileparts(here);
            tc.Examples = fullfile(tc.RepoRoot, 'examples');

            analyzerSrc = fullfile(tc.RepoRoot, 'ABL_CDI_Analyzer', 'src');
            patlogSrc   = fullfile(tc.RepoRoot, 'PatLogGUI', 'src');
            addpath(analyzerSrc);
            addpath(patlogSrc);
            tc.addTeardown(@() rmpath(analyzerSrc));
            tc.addTeardown(@() rmpath(patlogSrc));

            tc.AblFile        = fullfile(tc.Examples, 'PatLog_export.csv');
            tc.CdiFile2026027 = fullfile(tc.Examples, '2026027_cdi.log');
            tc.CdiFile2026007 = fullfile(tc.Examples, '2026007_cdi.log');

            tc.assertTrue(isfile(tc.AblFile),        'Missing examples/PatLog_export.csv');
            tc.assertTrue(isfile(tc.CdiFile2026027), 'Missing examples/2026027_cdi.log');
            tc.assertTrue(isfile(tc.CdiFile2026007), 'Missing examples/2026007_cdi.log');
        end
    end

    methods (Access = private)
        function app = newAnalyzer(tc)
            app = ABL_CDI_Analyzer;
            tc.addTeardown(@() delete(app));
        end

        function fig = launchPatLog(tc, file)
            evalc('PatLogGUI(file)');
            figs = findall(0, 'Type', 'figure');
            keep = arrayfun(@(f) contains(string(f.Name), ...
                'PatLog Blood Gas Analyzer'), figs);
            tc.assertTrue(any(keep), 'PatLogGUI window was not created');
            fig = figs(find(keep, 1));
            tc.addTeardown(@() delete(fig));
        end

        function pushButton(~, fig, label)
            btns = findall(fig, 'Type', 'uibutton');
            hit  = btns(arrayfun(@(b) contains(string(b.Text), label), btns));
            evalc('feval(hit(1).ButtonPushedFcn, hit(1), [])');
        end

        function txt = labelContaining(tc, fig, needle)
            labs = findall(fig, 'Type', 'uilabel');
            hit  = labs(arrayfun(@(l) contains(string(l.Text), needle), labs));
            tc.assertNotEmpty(hit, sprintf('No label containing "%s"', needle));
            txt = char(hit(1).Text);
        end

        function s = keptPairStats(tc, out)
            % Before/after agreement on the MAD-retained pairs, computed the
            % way Export Figures does for its Bland-Altman panels.
            x = out.stats.xABL; yRaw = out.stats.yCDI; yCor = out.model.yCorrected;
            tc.assertEqual(numel(yCor), numel(x), 'Corrected series and pairs differ in length');
            keep = ~isnan(yCor) & ~isnan(x);
            dB = yRaw(keep) - x(keep); dA = yCor(keep) - x(keep);
            s.n = nnz(keep);
            s.biasBefore = mean(dB); s.sdBefore = std(dB);
            s.biasAfter  = mean(dA); s.sdAfter  = std(dA);
            s.loaBefore = sprintf('[%.3f, %.3f]', s.biasBefore - 1.96*s.sdBefore, s.biasBefore + 1.96*s.sdBefore);
            s.loaAfter  = sprintf('[%.3f, %.3f]', s.biasAfter  - 1.96*s.sdAfter,  s.biasAfter  + 1.96*s.sdAfter);
        end
    end

    methods (Test)

        % ---------- WALKTHROUGH.md section 2: pH, dataset 2026027 ----------

        function section2_beforeCorrection(tc)
            % Steps 4-6, including the Fit Window / Auto step.
            app = tc.newAnalyzer();
            out = app.runWorkflow(tc.AblFile, tc.CdiFile2026027, '2026027', 'pH', ...
                'TimeTolerance', 5, 'FitWindowAuto', true);

            tc.verifyEqual(out.nPairsText, 'N Pairs: 67 (window) / 68 total');
            tc.verifyEqual(out.biasText,   'Bias: 0.047');
            tc.verifyEqual(out.sdText,     'SD: 0.095');
            tc.verifyEqual(out.loaText,    '95% LoA: [-0.140, 0.234]');
            tc.verifyEqual(out.rText,      'r = 0.1700');
            tc.verifyEqual(app.FitWindowStartEdit.Value, '14.04.2026 16:19');
            tc.verifyEqual(app.FitWindowEndEdit.Value,   '14.04.2026 21:03');
        end

        function section2_fitWindowIsRequired(tc)
            % Without the Fit Window step, the §2 values must not appear.
            app = tc.newAnalyzer();
            out = app.runWorkflow(tc.AblFile, tc.CdiFile2026027, '2026027', 'pH', ...
                'TimeTolerance', 5, 'FitWindowAuto', false, ...
                'CorrectionMethod', 'Auto (Best Model)');

            tc.verifyEqual(out.nPairsText, 'N Pairs: 68');
            tc.verifyNotEqual(out.biasText, 'Bias: 0.047');
            % ...and the values WALKTHROUGH.md says you get instead
            tc.verifyEqual(out.biasText, 'Bias: 0.049');
            tc.verifyEqual(out.sdText,   'SD: 0.096');
            tc.verifyEqual(out.loaText,  '95% LoA: [-0.139, 0.236]');
            tc.verifySubstring(out.formula, '1.8666*ABL -6.2502');
            tc.verifySubstring(out.formula, '68/68');
            tc.verifySubstring(out.formula, '62 robust');
            tc.verifySubstring(out.qualityText, 'SD ▼12.8%');
        end

        function section2_autoCorrection(tc)
            app = tc.newAnalyzer();
            out = app.runWorkflow(tc.AblFile, tc.CdiFile2026027, '2026027', 'pH', ...
                'TimeTolerance', 5, 'FitWindowAuto', true, ...
                'CorrectionMethod', 'Auto (Best Model)');

            tc.verifyEqual(out.model.autoWinner, 'Weighted Deming (Linnet, tuned λ)');
            tc.verifyEqual(out.model.slope,     1.7395, 'AbsTol', 5e-4);
            tc.verifyEqual(out.model.intercept, -5.3238, 'AbsTol', 5e-4);
            tc.verifyEqual(out.model.lam,        0.10,   'AbsTol', 1e-9);

            tc.verifySubstring(out.formula, '1.7395*ABL -5.3238');
            tc.verifySubstring(out.formula, 'Fitted on 67/68 window pairs, 61 robust');

            tc.verifyEqual(out.afterBiasText, 'After Correction: Bias=-0.0002');
            tc.verifyEqual(out.afterSDText,   'SD=0.0426 (on 61 kept pairs)');
            tc.verifyEqual(out.model.r_new,   0.4711, 'AbsTol', 5e-4);

            % Same 61 retained pairs before and after (Export Figures panels)
            s = tc.keptPairStats(out);
            tc.verifyEqual(s.n, 61);
            tc.verifyEqual(sprintf('%.4f', s.biasBefore), '0.0729');
            tc.verifyEqual(sprintf('%.4f', s.sdBefore),   '0.0475');
            tc.verifyEqual(s.loaBefore, '[-0.020, 0.166]');
            tc.verifyEqual(s.loaAfter,  '[-0.084, 0.083]');
            tc.verifySubstring(out.qualityText, 'BIAS + SD REDUCED');
            tc.verifySubstring(out.qualityText, 'SD ▼10.3%');

            % Reference line: mean of the training ABL values, ignoring the CDI
            tc.verifyEqual(sprintf('%.4f', out.model.autoReferenceRMSE), '0.0935');
            tc.verifyTrue(any(contains(out.model.autoRankText, 'RMSE=0.0935  (ignores the CDI)')));
            tc.verifyTrue(any(contains(out.model.autoRankText, 'No candidate beats the reference')));
        end

        function section2_tieBreakerEngages(tc)
            % The top two candidates fall inside the 1% RMSE band, so the
            % limits-of-agreement tie-breaker decides.
            app = tc.newAnalyzer();
            out = app.runWorkflow(tc.AblFile, tc.CdiFile2026027, '2026027', 'pH', ...
                'TimeTolerance', 5, 'FitWindowAuto', true, ...
                'CorrectionMethod', 'Auto (Best Model)');

            sorted = sort(out.model.autoRMSE);
            gapPct = 100 * (sorted(2) - sorted(1)) / sorted(1);
            tc.verifyLessThan(gapPct, 1, ...
                'Expected the top two candidates within 1% RMSE for pH/2026027');

            [~, order] = sort(out.model.autoRMSE);
            tc.verifyEqual(sprintf('%.4f', sorted(1)), '0.0953');
            tc.verifyEqual(sprintf('%.4f', sorted(2)), '0.0957');
            tc.verifyEqual(sprintf('%.2f', gapPct),    '0.42');
            tc.verifyEqual(sprintf('%.4f', out.model.autoLoASpan(order(1))), '0.3628');
            tc.verifyEqual(sprintf('%.4f', out.model.autoLoASpan(order(2))), '0.3646');
        end

        % ---------- WALKTHROUGH.md section 3: pO2, dataset 2026007 ----------

        function section3_beforeCorrection(tc)
            app = tc.newAnalyzer();
            out = app.runWorkflow(tc.AblFile, tc.CdiFile2026007, '2026007', 'pO2', ...
                'TimeTolerance', 5);

            tc.verifyEqual(out.nPairsText, 'N Pairs: 10');
            tc.verifyEqual(out.biasText,   'Bias: 48.460');
            tc.verifyEqual(out.sdText,     'SD: 277.922');
            tc.verifyEqual(out.loaText,    '95% LoA: [-496.267, 593.187]');
            tc.verifyEqual(out.rText,      'r = -0.1088');
        end

        function section3_autoCorrection(tc)
            app = tc.newAnalyzer();
            out = app.runWorkflow(tc.AblFile, tc.CdiFile2026007, '2026007', 'pO2', ...
                'TimeTolerance', 5, 'CorrectionMethod', 'Auto (Best Model)');

            tc.verifyEqual(out.model.autoWinner, 'Hybrid (Time-Series + Deming)');
            tc.verifyEqual(min(out.model.autoRMSE), 205.9690, 'AbsTol', 5e-4);
            tc.verifyEqual(out.model.slope,     -2.2513,   'AbsTol', 5e-4);
            tc.verifyEqual(out.model.intercept,  698.8642, 'AbsTol', 5e-3);
            tc.verifyEqual(out.model.lam,        0.50,     'AbsTol', 1e-9);
            tc.verifyEqual(out.model.tau_rise,   8.0,      'AbsTol', 1e-9);
            tc.verifyEqual(out.model.tau_fall,   8.0,      'AbsTol', 1e-9);
            tc.verifyEqual(out.model.w1,         32,       'AbsTol', 1e-9);

            tc.verifySubstring(out.formula, '(CDI_fast - 698.8642) / -2.2513');
            tc.verifySubstring(out.formula, 'Fitted on 10/10 window pairs, 10 robust');

            tc.verifyEqual(out.afterBiasText, 'After Correction: Bias=-81.9161');
            tc.verifyEqual(out.afterSDText,   'SD=180.2320 (on 10 kept pairs)');
            tc.verifyEqual(out.model.r_new,   0.1573, 'AbsTol', 5e-4);
            tc.verifySubstring(out.qualityText, 'SD ▼35.2%');
            % the bias grows in magnitude, so the label is SD REDUCED only
            tc.verifyTrue(startsWith(out.qualityText, 'SD REDUCED'));

            s = tc.keptPairStats(out);
            tc.verifyEqual(s.n, 10);
            tc.verifyEqual(s.loaAfter, '[-435.171, 271.339]');

            % Reference line: mean of the training ABL values, ignoring the CDI
            tc.verifyEqual(sprintf('%.4f', out.model.autoReferenceRMSE), '173.8841');
            tc.verifyTrue(any(contains(out.model.autoRankText, 'RMSE=173.8841  (ignores the CDI)')));
            tc.verifyTrue(any(contains(out.model.autoRankText, 'No candidate beats the reference')));
        end

        function section3_noTieBreaker(tc)
            % Unlike the pH dataset, the top two candidates are more than 1%
            % apart, so the tie-breaker must not engage.
            app = tc.newAnalyzer();
            out = app.runWorkflow(tc.AblFile, tc.CdiFile2026007, '2026007', 'pO2', ...
                'TimeTolerance', 5, 'CorrectionMethod', 'Auto (Best Model)');

            sorted = sort(out.model.autoRMSE);
            gapPct = 100 * (sorted(2) - sorted(1)) / sorted(1);
            tc.verifyGreaterThan(gapPct, 1, ...
                'Expected the top two candidates further than 1% apart for pO2/2026007');
            tc.verifyEqual(sprintf('%.4f', sorted(2)), '213.2436');
            tc.verifyEqual(sprintf('%.1f', gapPct),    '3.5');
        end

        function section3_allPatientsGivesSameResult(tc)
            % The walkthrough's note: this CDI recording only overlaps patient
            % 2026007, so leaving the dropdown on "All Patients" is equivalent.
            appA = tc.newAnalyzer();
            selected = appA.runWorkflow(tc.AblFile, tc.CdiFile2026007, '2026007', 'pO2', ...
                'TimeTolerance', 5);
            appB = tc.newAnalyzer();
            allPat = appB.runWorkflow(tc.AblFile, tc.CdiFile2026007, 'All Patients', 'pO2', ...
                'TimeTolerance', 5);

            tc.verifyEqual(allPat.nPairsText, selected.nPairsText);
            tc.verifyEqual(allPat.stats.bias, selected.stats.bias, 'AbsTol', 1e-12);
            tc.verifyEqual(allPat.stats.sd,   selected.stats.sd,   'AbsTol', 1e-12);
        end

        function looCvTuningIsFoldIndependent(tc)
            % Hybrid's tuning statistics - the derivative clip and the roughness
            % reference - are computed per training fold. Dropping an interior
            % pair leaves the fold's time span unchanged, so leave-one-out needs
            % only three distinct spans; this checks the winner and its RMSE are
            % stable across repeated runs, which they would not be if a
            % fold-dependent quantity were being read from the wrong span.
            first = [];
            for k = 1:2
                app = tc.newAnalyzer();
                out = app.runWorkflow(tc.AblFile, tc.CdiFile2026007, '2026007', 'pO2', ...
                    'TimeTolerance', 5, 'CorrectionMethod', 'Auto (Best Model)');
                if isempty(first)
                    first = out.model.autoRMSE;
                    tc.verifyEqual(out.model.autoWinner, 'Hybrid (Time-Series + Deming)');
                else
                    tc.verifyEqual(out.model.autoRMSE, first, 'AbsTol', 1e-12, ...
                        'LOO-CV RMSE must be reproducible run to run');
                end
                tc.verifyFalse(any(isnan(out.model.autoRMSE)), ...
                    'every candidate must score - a NaN means folds are erroring out');
            end
        end

        % ---------- WALKTHROUGH.md section 4: batch analysis ----------

        function section4_batchReproducesExamples(tc)
            % runBatch on examples/batch_example.csv gives the section 2 and 3
            % results in one summary table, and writes it to disk.
            outDir = tempname;
            mkdir(outDir);
            tc.addTeardown(@() rmdir(outDir, 's'));
            outFile = fullfile(outDir, 'summary.csv');
            evalc('S = ABL_CDI_Analyzer.runBatch(fullfile(tc.Examples, ''batch_example.csv''), outFile);');

            tc.verifyEqual(height(S), 2);
            tc.verifyTrue(all(S.Status == "ok"));
            tc.verifyEqual(height(readtable(outFile)), 2);

            pH = S(S.Parameter == "pH", :);
            tc.verifyEqual([pH.NPairs pH.NTotal pH.NKept], [67 68 61]);
            tc.verifyEqual(pH.Model, "Weighted Deming (Linnet, tuned λ)");
            tc.verifyEqual(sprintf('%.4f %.4f %.4f', pH.LOOCV_RMSE, pH.RunnerUp_RMSE, pH.Reference_RMSE), ...
                '0.0953 0.0957 0.0935');
            tc.verifyEqual(pH.BeatsReference, "no");
            tc.verifyEqual(sprintf('%.4f %.4f', pH.BiasAfter, pH.SDAfter), '-0.0002 0.0426');
            tc.verifyEqual(sprintf('[%.3f, %.3f]', pH.LoALowAfter, pH.LoAHighAfter), '[-0.084, 0.083]');
            tc.verifyEqual(sprintf('%.1f', pH.SDReductionPct), '10.3');
            tc.verifyEqual(pH.Label, "BIAS + SD REDUCED");

            pO2 = S(S.Parameter == "pO2", :);
            tc.verifyEqual([pO2.NPairs pO2.NTotal pO2.NKept], [10 10 10]);
            tc.verifyEqual(pO2.Model, "Hybrid (Time-Series + Deming)");
            tc.verifyEqual(sprintf('%.4f %.4f %.4f', pO2.LOOCV_RMSE, pO2.RunnerUp_RMSE, pO2.Reference_RMSE), ...
                '205.9690 213.2436 173.8841');
            tc.verifyEqual(pO2.BeatsReference, "no");
            tc.verifyEqual(sprintf('%.4f %.4f', pO2.BiasAfter, pO2.SDAfter), '-81.9161 180.2320');
            tc.verifyEqual(sprintf('[%.3f, %.3f]', pO2.LoALowAfter, pO2.LoAHighAfter), '[-435.171, 271.339]');
            tc.verifyEqual(sprintf('%.1f', pO2.SDReductionPct), '35.2');
            tc.verifyEqual(pO2.Label, "SD REDUCED");
        end

        function section4_batchContinuesAfterFailure(tc)
            % A recording that cannot be read is reported, not fatal.
            bad = table("missing.csv", "2026007_cdi.log", "2026007", "pO2", ...
                'VariableNames', {'ABLFile', 'CDIFile', 'PatientID', 'Parameter'});
            evalc('S = ABL_CDI_Analyzer.runBatch(bad);');
            tc.verifyEqual(height(S), 1);
            tc.verifyTrue(startsWith(S.Status, "error:"));
        end

        % ---------- feature regressions ----------

        function exportResultsWritesFiles(tc)
            % Export Results writes scalar statistics and paired data separately.
            app = tc.newAnalyzer();
            app.runWorkflow(tc.AblFile, tc.CdiFile2026027, '2026027', 'pH', ...
                'TimeTolerance', 5, 'FitWindowAuto', true);

            outDir = tempname;
            mkdir(outDir);
            tc.addTeardown(@() rmdir(outDir, 's'));

            app.exportResultsTo(fullfile(outDir, 'results.xlsx'));
            tc.verifyTrue(isfile(fullfile(outDir, 'results.xlsx')), ...
                'XLSX export produced no file');

            app.exportResultsTo(fullfile(outDir, 'results.csv'));
            tc.verifyTrue(isfile(fullfile(outDir, 'results_statistics.csv')), ...
                'CSV export produced no statistics file');
            tc.verifyTrue(isfile(fullfile(outDir, 'results_paireddata.csv')), ...
                'CSV export produced no paired-data file');
        end

        function patLogKeepsEighteenEssentialColumns(tc)
            fig = tc.launchPatLog(tc.AblFile);
            tc.pushButton(fig, 'Process Clinical Data');
            info = tc.labelContaining(fig, 'essential columns');
            tc.verifySubstring(info, 'Kept 18 essential columns');
            tc.verifySubstring(info, 'Removed 218 unnecessary columns');
            % cleaning selects columns only: every data row is still shown
            tc.verifySubstring(tc.labelContaining(fig, 'Displaying'), 'Displaying 1907 rows');
        end

        function patLogReadsLatin1Encoding(tc)
            % examples/PatLog_export.csv is ISO-8859-1: the degree sign in the
            % temperature header is a single 0xB0 byte, and the column name must
            % read "T (°C)".
            fig = tc.launchPatLog(tc.AblFile);
            tc.pushButton(fig, 'Process Clinical Data');
            info = tc.labelContaining(fig, 'essential columns');

            tc.verifyFalse(contains(info, char(65533)), ...
                'Temperature column name contains a Unicode replacement character');
            tc.verifySubstring(info, ['T (' char(176) 'C)']);
        end

        function toolboxDependenciesMatchDocumentation(tc)
            % README, CODE_METADATA.md and the per-tool READMEs state that the
            % Analyzer is core-MATLAB-only and PatLogGUI additionally needs the
            % Statistics and Machine Learning Toolbox.
            analyzer = fullfile(tc.RepoRoot, 'ABL_CDI_Analyzer', 'src', 'ABL_CDI_Analyzer.m');
            patlog   = fullfile(tc.RepoRoot, 'PatLogGUI', 'src', 'PatLogGUI.m');

            [~, prodA] = matlab.codetools.requiredFilesAndProducts(analyzer);
            tc.verifyEqual(sort(string({prodA.Name})), "MATLAB", ...
                'ABL_CDI_Analyzer is documented as core MATLAB only');

            [~, prodP] = matlab.codetools.requiredFilesAndProducts(patlog);
            tc.verifyEqual(sort(string({prodP.Name})), ...
                sort(["MATLAB", "Statistics and Machine Learning Toolbox"]), ...
                'PatLogGUI dependencies differ from what the READMEs state');
        end
    end
end
