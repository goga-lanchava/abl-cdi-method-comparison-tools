classdef ComponentTest < matlab.unittest.TestCase
    % Deterministic component tests for ABL-CDI Analyzer: parsing, input
    % validation, temporal pairing, MAD filtering and one numerical
    % correction with an independently computed expected value.
    %
    % Run with:   runtests('tests')

    properties
        RepoRoot
        Examples
        AblFile
        CdiFile2026027
        TempDir
    end

    methods (TestClassSetup)
        function setupPaths(tc)
            here = fileparts(which('ComponentTest'));
            tc.RepoRoot = fileparts(here);
            tc.Examples = fullfile(tc.RepoRoot, 'examples');

            src = fullfile(tc.RepoRoot, 'ABL_CDI_Analyzer', 'src');
            addpath(src);
            tc.addTeardown(@() rmpath(src));

            tc.AblFile        = fullfile(tc.Examples, 'PatLog_export.csv');
            tc.CdiFile2026027 = fullfile(tc.Examples, '2026027_cdi.log');

            tc.TempDir = tempname;
            mkdir(tc.TempDir);
            tc.addTeardown(@() rmdir(tc.TempDir, 's'));
        end
    end

    methods (Access = private)
        function app = newApp(tc)
            app = ABL_CDI_Analyzer;
            tc.addTeardown(@() delete(app));
        end

        function f = writeFile(tc, name, lines)
            f = fullfile(tc.TempDir, name);
            fid = fopen(f, 'w');
            tc.assertNotEqual(fid, -1);
            fprintf(fid, '%s\n', lines{:});
            fclose(fid);
        end
    end

    methods (Test)

        % ---------------- parsing a known ABL export ----------------

        function ablParserReadsKnownFile(tc)
            app = tc.newApp();
            [tbl, ids] = app.readABL(tc.AblFile);

            tc.verifyEqual(height(tbl), 1907);
            tc.verifyTrue(all(ismember({'Time','PatientID','pH','pO2','pCO2'}, ...
                tbl.Properties.VariableNames)));
            tc.verifyTrue(ismember('2026027', ids));
            tc.verifyTrue(ismember('2026007', ids));
            tc.verifyClass(tbl.Time, 'datetime');
            tc.verifyTrue(all(~isnat(tbl.Time)));
        end

        % ---------------- parsing a known CDI log ----------------

        function cdiParserReadsKnownFile(tc)
            app = tc.newApp();
            tbl = app.readCDI(tc.CdiFile2026027);

            tc.verifyEqual(height(tbl), 3235);
            tc.verifyEqual(tbl.Properties.VariableNames, ...
                {'Time','pH','pCO2','pO2','TEMP','HCO3','BE','sO2','K+'});
            tc.verifyEqual(min(tbl.Time), datetime(2026,4,14,16,16,40));
            tc.verifyEqual(max(tbl.Time), datetime(2026,4,14,21,38,47));
        end

        function cdiParserDropsEntirelyMissingColumns(tc)
            % The log carries 18 fixed fields; only those with at least one
            % non-placeholder value survive. This recording has 8 parameters.
            app = tc.newApp();
            tbl = app.readCDI(tc.CdiFile2026027);
            tc.verifyEqual(width(tbl) - 1, 8);
            tc.verifyFalse(ismember('VO2', tbl.Properties.VariableNames));
        end

        function cdiParserSkipsMalformedLines(tc)
            % Lines without a bracketed date, without a device time, or with
            % only placeholders must be discarded; valid lines must survive.
            f = tc.writeFile('synthetic_cdi.log', { ...
                'garbage header line with no bracket', ...
                '[2026-04-14 16:00:00.000] 16:00:00,7.400,40.0,100.0,37.0,24.0,0.0,98.0,4.0,,,,,,,,,', ...
                '[2026-04-14 16:01:00.000] notatime,7.410,40.1,101.0,37.0,24.0,0.0,98.1,4.0,,,,,,,,,', ...
                '[2026-04-14 16:02:00.000] 16:02:00,%,---,--,-.-,,,,,,,,,,,,,', ...
                '[2026-04-14 16:03:00.000] 16:03:00,7.420,40.2,102.0,37.0,24.0,0.0,98.2,4.0,,,,,,,,,'});
            app = tc.newApp();
            tbl = app.readCDI(f);

            tc.verifyEqual(height(tbl), 2, ...
                'Only the two fully valid rows should be kept');
            tc.verifyEqual(tbl.pH, [7.400; 7.420], 'AbsTol', 1e-9);
        end

        % ---------------- input validation ----------------

        function ablParserRejectsMalformedTimestamps(tc)
            % The first column must parse as d.M.yyyy HH:mm.
            f = tc.writeFile('bad_time.csv', { ...
                'Time,Patient Id,pH', ...
                'not-a-timestamp,2026027,7.40', ...
                'also-not-a-time,2026027,7.41'});
            app = tc.newApp();
            tc.verifyError(@() app.readABL(f), ?MException);
        end

        function ablParserToleratesMissingPatientColumn(tc)
            % Documents actual behaviour: a missing patient-ID column is not
            % fatal - the file parses and the patient list comes back empty.
            f = tc.writeFile('no_patient.csv', { ...
                'Time,pH,pO2', ...
                '14.04.2026 16:20,7.40,100', ...
                '14.04.2026 16:30,7.41,101'});
            app = tc.newApp();
            [tbl, ids] = app.readABL(f);

            tc.verifyEqual(height(tbl), 2);
            tc.verifyEmpty(ids);
        end

        % ---------------- temporal pairing ----------------

        function temporalPairingRespectsTolerance(tc)
            % Tightening the tolerance from 5 to 1 minute drops exactly the one
            % ABL draw whose nearest CDI sample is more than a minute away.
            appWide = tc.newApp();
            wide = appWide.runWorkflow(tc.AblFile, tc.CdiFile2026027, ...
                '2026027', 'pH', 'TimeTolerance', 5);
            appTight = tc.newApp();
            tight = appTight.runWorkflow(tc.AblFile, tc.CdiFile2026027, ...
                '2026027', 'pH', 'TimeTolerance', 1);

            tc.verifyEqual(wide.stats.Ntotal, 68);
            tc.verifyEqual(tight.stats.Ntotal, 67);
        end

        % ---------------- MAD filtering ----------------

        function madFilterExcludesGrossOutlierOnly(tc)
            % Differences 1..10 plus one gross outlier. Median difference is
            % 5.5, MAD is 2.5, so the retention band is 5.5 +/- 11.25 and only
            % the outlier falls outside it.
            app = tc.newApp();
            x = zeros(11, 1);
            y = [(1:10)'; 1000];
            mask = app.madRetainMask(x, y);

            tc.verifyEqual(mask, [true(10,1); false]);
        end

        function madFilterKeepsEverythingWhenMadIsZero(tc)
            % Degenerate case: identical differences give MAD = 0, and the
            % filter must not then reject every pair.
            app = tc.newApp();
            mask = app.madRetainMask(zeros(6,1), 3*ones(6,1));
            tc.verifyEqual(mask, true(6,1));
        end

        % ---------------- a numerical correction with a known answer ----------------

        function biasCorrectionMatchesHandComputedOffset(tc)
            % Bias Correction subtracts the mean paired difference computed on
            % the MAD-retained pairs. Recompute that independently and check
            % both the fitted offset and the corrected series against it.
            app = tc.newApp();
            out = app.runWorkflow(tc.AblFile, tc.CdiFile2026027, '2026027', 'pH', ...
                'TimeTolerance', 5, 'FitWindowAuto', true, ...
                'CorrectionMethod', 'Bias Correction');

            x = out.stats.xABL;
            y = out.stats.yCDI;
            keep = app.madRetainMask(x, y);
            expectedOffset = mean(y(keep) - x(keep), 'omitnan');

            tc.verifyEqual(out.model.bias, expectedOffset, 'AbsTol', 1e-10);

            corrected = out.model.yCorrected;
            tc.verifyEqual(corrected(keep), y(keep) - expectedOffset, 'AbsTol', 1e-10);
            tc.verifyTrue(all(isnan(corrected(~keep))), ...
                'MAD-excluded pairs must be NaN in the corrected series');

            residualBias = mean(corrected(keep) - x(keep), 'omitnan');
            tc.verifyEqual(residualBias, 0, 'AbsTol', 1e-10);
        end
    end
end
