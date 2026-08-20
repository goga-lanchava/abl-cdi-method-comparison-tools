classdef ABL_CDI_Analyzer < matlab.apps.AppBase

    properties (Access = public)
        UIFigure             matlab.ui.Figure
        GridLayout           matlab.ui.container.GridLayout
        LeftPanel            matlab.ui.container.Panel
        LeftGrid             matlab.ui.container.GridLayout

        LoadABLButton        matlab.ui.control.Button
        LoadCDIButton        matlab.ui.control.Button
        PatientIDLabel       matlab.ui.control.Label
        PatientIDDropDown    matlab.ui.control.DropDown
        ParamDropDown        matlab.ui.control.DropDown
        ParamLabel           matlab.ui.control.Label
        TimeToleranceSpinner matlab.ui.control.Spinner
        TimeToleranceLabel   matlab.ui.control.Label
        
        TimeShiftSpinner     matlab.ui.control.Spinner  
        TimeShiftLabel       matlab.ui.control.Label    
        AutoShiftButton      matlab.ui.control.Button   
        
        AnalyzeButton        matlab.ui.control.Button

        CorrectionPanel      matlab.ui.container.Panel
        CorrectionGrid       matlab.ui.container.GridLayout
        CorrectionMethodDropDown matlab.ui.control.DropDown
        CorrectionMethodLabel    matlab.ui.control.Label
        DemingLambdaLabel        matlab.ui.control.Label
        DemingLambdaEditField    matlab.ui.control.NumericEditField
        TauLabel                 matlab.ui.control.Label
        TauSpinner               matlab.ui.control.Spinner
        TauFallLabel             matlab.ui.control.Label
        TauFallSpinner           matlab.ui.control.Spinner
        SmoothW1Label            matlab.ui.control.Label
        SmoothW1Spinner          matlab.ui.control.Spinner
        AutoTuneButton           matlab.ui.control.Button 
        ApplyCorrectionButton    matlab.ui.control.Button
        ExportCorrectedButton    matlab.ui.control.Button
        ShowFormulaButton        matlab.ui.control.Button
        ComparePlotsButton       matlab.ui.control.Button 
        CorrectionStatusLabel    matlab.ui.control.Label
        FormulaLabel             matlab.ui.control.Label
        FormulaTextArea          matlab.ui.control.TextArea
        SmallNWarningLabel       matlab.ui.control.Label

        ExportButton         matlab.ui.control.Button
        ExportFigureButton   matlab.ui.control.Button
        StatusLabel          matlab.ui.control.Label

        StatsPanel           matlab.ui.container.Panel
        StatsGrid            matlab.ui.container.GridLayout
        NPairsLabel          matlab.ui.control.Label
        BiasLabel            matlab.ui.control.Label
        SDLabel              matlab.ui.control.Label
        LOALabel             matlab.ui.control.Label
        CorrelationLabel     matlab.ui.control.Label
        RegressionLabel      matlab.ui.control.Label
        ImprovedBiasLabel    matlab.ui.control.Label
        ImprovedSDLabel      matlab.ui.control.Label
        SlopeWarningLabel    matlab.ui.control.Label
        CorrQualityLabel     matlab.ui.control.Label

        RightPanel           matlab.ui.container.Panel
        RightGrid            matlab.ui.container.GridLayout
        TimeAxes             matlab.ui.control.UIAxes
        CorrelationAxes      matlab.ui.control.UIAxes
        BlandAltmanAxes      matlab.ui.control.UIAxes
        CorrectionAxes       matlab.ui.control.UIAxes
        ShowCorrectedSwitch  matlab.ui.control.CheckBox
        ShowOriginalCDISwitch matlab.ui.control.CheckBox
        KeepCorrectionSwitch matlab.ui.control.CheckBox
        CorrViewSwitch       matlab.ui.control.CheckBox  
        BAViewSwitch         matlab.ui.control.CheckBox  
        CorrectedTrendLine   % handle to plotted corrected line
        OriginalCDILine      % handle to original CDI line on time plot
        KeptCorrectionLines  % cell array of kept correction line handles
        PhysioGateMarkers    % (reserved, unused)

        % --- Fitting window controls ---
        FitWindowPanel       matlab.ui.container.Panel
        FitWindowCheckBox    matlab.ui.control.CheckBox
        FitWindowStartLabel  matlab.ui.control.Label
        FitWindowStartEdit   matlab.ui.control.EditField
        FitWindowEndLabel    matlab.ui.control.Label
        FitWindowEndEdit     matlab.ui.control.EditField
        FitWindowAutoButton  matlab.ui.control.Button
        StabilityScoreLabel  matlab.ui.control.Label
        FitWindowShade       % patch handle - shaded region on time plot
    end

    properties (Access = private)
        ABL_Table            table = table()
        ABL_Table_Full       table = table()   
        ABL_PatientIDs       cell = {}         
        CDI_Table            table = table()
        CDI_Corrected_Table  table = table()
        Aligned_Data         timetable = timetable()
        CurrentParam         char = ''
        Stats                struct = struct()
        CorrectionModel      struct = struct()
        KeptCorrectionData   cell = {}  
    end

    properties (Access = private, Constant)
        RefRanges = struct()
    end

    methods (Access = private)

        % --- MAIN FILTERING FUNCTION ---
        % Uses Median Absolute Deviation (MAD) to filter extreme outliers.
        function cleanMask = robustCleanMask(~, xABL, yCDI)
            diffs = yCDI - xABL;
            med   = median(diffs, 'omitnan');
            mad_  = median(abs(diffs - med), 'omitnan');
            if mad_ > 0
                cleanMask = abs(diffs - med) <= (4.5 * mad_);
            else
                cleanMask = true(size(diffs));
            end
        end

        % --- WEIGHTED DEMING REGRESSION (LINNET ALGORITHM) ---
        % Iteratively computes weighted Deming regression where weight w_i = 1 / u_hat_i^2.
        function [slope, intercept] = fitWeightedDeming(~, x, y, lambda)
            if nargin < 4 || isempty(lambda), lambda = 1.0; end
            n = numel(x);
            if n < 2
                slope = 1.0; intercept = 0.0; return;
            end

            xm = mean(x, 'omitnan'); ym = mean(y, 'omitnan');
            sxx = sum((x - xm).^2, 'omitnan') / max(n - 1, 1);
            syy = sum((y - ym).^2, 'omitnan') / max(n - 1, 1);
            sxy = sum((x - xm) .* (y - ym), 'omitnan') / max(n - 1, 1);
            
            if abs(2 * sxy) > 1e-10
                slope = (syy - lambda * sxx + sqrt((syy - lambda * sxx)^2 + 4 * lambda * sxy^2)) / (2 * sxy);
            else
                slope = 1.0;
            end
            if isnan(slope) || isinf(slope) || abs(slope) < 1e-4, slope = 1.0; end
            intercept = ym - slope * xm;

            for iter = 1:50
                prev_slope = slope;
                prev_intercept = intercept;

                u_hat = 0.5 * (x + (y - intercept) / max(slope, 1e-5));
                u_hat(abs(u_hat) < 1e-4) = 1e-4;
                w = 1 ./ (u_hat.^2);
                w = w / sum(w); 

                x_w = sum(w .* x);
                y_w = sum(w .* y);

                dx = x - x_w;
                dy = y - y_w;

                sxx_w = sum(w .* dx.^2);
                syy_w = sum(w .* dy.^2);
                sxy_w = sum(w .* dx .* dy);

                if abs(2 * sxy_w) > 1e-10
                    slope = (syy_w - lambda * sxx_w + sqrt((syy_w - lambda * sxx_w)^2 + 4 * lambda * sxy_w^2)) / (2 * sxy_w);
                else
                    slope = 1.0;
                end
                if isnan(slope) || isinf(slope) || abs(slope) < 1e-4, slope = 1.0; end
                intercept = y_w - slope * x_w;

                if abs(slope - prev_slope) < 1e-5 && abs(intercept - prev_intercept) < 1e-5
                    break;
                end
            end
        end

        % --- CAUSAL DYNAMIC RESPONSE FILTER ---
        function cdi_fast = computeAsymmetricFastCDI(~, fullCDIVals, fullCDITime, w1, tau_rise, tau_fall)
            w_sm = max(0, w1 - 1);
            smoothed = movmean(fullCDIVals, [w_sm 0], 'omitnan');
            if tau_rise == 0 && tau_fall == 0
                cdi_fast = smoothed;
            else
                dt = minutes(diff(fullCDITime));
                dt(dt <= 0 | ~isfinite(dt)) = 0.01;
                dy = diff(smoothed);
                if isempty(dy)
                    deriv = zeros(size(fullCDIVals));
                else
                    deriv = [0; dy ./ dt];
                    deriv = movmean(deriv, [w_sm 0], 'omitnan');
                end
                validD = deriv(~isnan(deriv) & ~isinf(deriv));
                if ~isempty(validD)
                    dlim = prctile(abs(validD), 95);
                    if dlim == 0 || isnan(dlim), dlim = 10; end
                else
                    dlim = 10;
                end
                deriv = max(min(deriv, dlim), -dlim);
                
                tau_vec = zeros(size(deriv));
                tau_vec(deriv >= 0) = tau_rise;
                tau_vec(deriv < 0)  = tau_fall;
                
                cdi_fast = smoothed + tau_vec .* deriv;
            end
            cdi_fast(isnan(cdi_fast)) = fullCDIVals(isnan(cdi_fast));
            cdi_fast = movmean(cdi_fast, [2 0], 'omitnan');
        end

        function [cvPct, nPairsInWindow] = computeStabilityScore(app, fitStart, fitEnd, param, pairedTimes)
            cdiTimes = app.CDI_Table.Time + minutes(app.TimeShiftSpinner.Value);
            if ismember(param, app.CDI_Table.Properties.VariableNames)
                cdiVals = app.CDI_Table.(param);
            else
                cvPct = NaN; nPairsInWindow = 0; return;
            end
            winMask = cdiTimes >= fitStart & cdiTimes <= fitEnd & ~isnan(cdiVals);
            winVals = cdiVals(winMask);
            if numel(winVals) < 3 || mean(abs(winVals),'omitnan') < 1e-6
                cvPct = NaN;
            else
                cvPct = 100 * std(winVals,'omitnan') / mean(abs(winVals),'omitnan');
            end
            nPairsInWindow = sum(pairedTimes >= fitStart & pairedTimes <= fitEnd);
        end

        function winMask = getFitWindowMask(app, pairedTimes)
            if ~app.FitWindowCheckBox.Value
                winMask = true(size(pairedTimes));
                return;
            end
            try
                fitStart = datetime(app.FitWindowStartEdit.Value, ...
                    'InputFormat', 'dd.MM.yyyy HH:mm');
                fitEnd   = datetime(app.FitWindowEndEdit.Value, ...
                    'InputFormat', 'dd.MM.yyyy HH:mm');
                winMask = pairedTimes >= fitStart & pairedTimes <= fitEnd;
                if sum(winMask) < 2
                    winMask = true(size(pairedTimes));
                    uialert(app.UIFigure, ...
                        'Fitting window contains fewer than 2 pairs - using all pairs instead.', ...
                        'Window Too Narrow');
                end
            catch
                winMask = true(size(pairedTimes));
            end
        end

        function drawFitWindowShade(app)
            if ~isempty(app.FitWindowShade) && isvalid(app.FitWindowShade)
                delete(app.FitWindowShade);
                app.FitWindowShade = [];
            end
            if ~app.FitWindowCheckBox.Value
                return;
            end
            try
                fitStart = datetime(app.FitWindowStartEdit.Value, ...
                    'InputFormat', 'dd.MM.yyyy HH:mm');
                fitEnd   = datetime(app.FitWindowEndEdit.Value, ...
                    'InputFormat', 'dd.MM.yyyy HH:mm');
            catch
                return;
            end
            yl = ylim(app.TimeAxes);
            hold(app.TimeAxes, 'on');
            app.FitWindowShade = patch(app.TimeAxes, ...
                [fitStart fitEnd fitEnd fitStart], ...
                [yl(1) yl(1) yl(2) yl(2)], ...
                [0.2 0.6 1.0], ...
                'FaceAlpha', 0.10, ...
                'EdgeColor', [0.2 0.6 1.0], ...
                'EdgeAlpha', 0.5, ...
                'LineStyle', '--', ...
                'DisplayName', 'Fitting window');
            uistack(app.FitWindowShade, 'bottom');
            hold(app.TimeAxes, 'off');
            legend(app.TimeAxes, 'Location', 'best');
        end

        function LoadABLButtonPushed(app, ~)
            [file, path] = uigetfile({'*.csv';'*.xlsx';'*.xls'}, 'Select ABL Data');
            if isequal(file,0), return; end
            app.StatusLabel.Text = 'Loading ABL...'; drawnow;
            try
                fullpath = fullfile(path, file);
                [tbl, patientIDs] = parseABL(app, fullpath);

                app.ABL_Table_Full = tbl;
                app.ABL_PatientIDs = patientIDs;
                app.ABL_Table = tbl;  
                
                app.PatientIDDropDown.Items = [{'All Patients'}, patientIDs(:)'];
                app.PatientIDDropDown.Value = 'All Patients';
                
                app.StatusLabel.Text = sprintf('ABL Loaded: %d samples, %d parameters', height(app.ABL_Table), width(app.ABL_Table)-2);

                if ~isempty(app.CDI_Table)
                    updateCommonParameters(app);
                else
                    cols = app.ABL_Table.Properties.VariableNames;
                    paramCols = setdiff(cols, {'Time', 'PatientID'});
                    if isempty(paramCols)
                        app.ParamDropDown.Items = {'No numeric parameter found'};
                        app.ParamDropDown.Value = app.ParamDropDown.Items{1};
                    else
                        app.ParamDropDown.Items = sort(paramCols);
                        app.ParamDropDown.Value = paramCols{1};
                    end
                end
            catch ME
                uialert(app.UIFigure, ['Error loading ABL: ' ME.message], 'Load Error');
                app.StatusLabel.Text = 'ABL Load Failed.';
            end
        end
        
        function PatientIDDropDownValueChanged(app, ~)
            selectedID = app.PatientIDDropDown.Value;
            
            if strcmp(selectedID, 'All Patients')
                app.ABL_Table = app.ABL_Table_Full;
            else
                if ismember('PatientID', app.ABL_Table_Full.Properties.VariableNames)
                    mask = strcmp(app.ABL_Table_Full.PatientID, selectedID);
                    app.ABL_Table = app.ABL_Table_Full(mask, :);
                else
                    app.ABL_Table = app.ABL_Table_Full;
                end
            end
            
            app.StatusLabel.Text = sprintf('ABL Filtered: %d samples', height(app.ABL_Table));
            
            if ~isempty(app.CDI_Table)
                updateCommonParameters(app);
            else
                cols = app.ABL_Table.Properties.VariableNames;
                paramCols = setdiff(cols, {'Time', 'PatientID'});
                if ~isempty(paramCols)
                    currentParam = app.ParamDropDown.Value;
                    app.ParamDropDown.Items = sort(paramCols);
                    if ismember(currentParam, paramCols)
                        app.ParamDropDown.Value = currentParam;
                    else
                        app.ParamDropDown.Value = paramCols{1};
                    end
                end
            end
        end

        function LoadCDIButtonPushed(app, ~)
            [file, path] = uigetfile({'*.csv';'*.txt';'*.log';'*.*'}, 'Select CDI Data');
            if isequal(file,0), return; end
            app.StatusLabel.Text = 'Loading CDI...'; drawnow;
            try
                fullpath = fullfile(path, file);
                tbl = parseCDI(app, fullpath);

                app.CDI_Table = tbl;
                app.StatusLabel.Text = sprintf('CDI Loaded: %d samples, %d parameters', height(app.CDI_Table), width(app.CDI_Table)-1);

                updateCommonParameters(app);
            catch ME
                uialert(app.UIFigure, ['Error loading CDI: ' ME.message], 'Load Error');
                app.StatusLabel.Text = 'CDI Load Failed.';
            end
        end
        
        function updateCommonParameters(app)
            if isempty(app.ABL_Table) || isempty(app.CDI_Table)
                return;
            end
            
            ablCols = setdiff(app.ABL_Table.Properties.VariableNames, {'PatientID', 'Time'});
            cdiCols = setdiff(app.CDI_Table.Properties.VariableNames, {'Time'});
            
            common = intersect(ablCols, cdiCols, 'stable');
            
            for i = 1:numel(ablCols)
                ablName = ablCols{i};
                ablClean = lower(regexprep(ablName, '[^a-zA-Z0-9]', ''));
                
                for j = 1:numel(cdiCols)
                    cdiName = cdiCols{j};
                    cdiClean = lower(regexprep(cdiName, '[^a-zA-Z0-9]', ''));
                    
                    if strcmp(ablClean, cdiClean) && ~ismember(ablName, common)
                        common{end+1} = ablName; 
                    end
                end
            end
            
            if ~isempty(common)
                common = sort(common);
                currentParam = app.ParamDropDown.Value;
                app.ParamDropDown.Items = common;
                if ismember(currentParam, common)
                    app.ParamDropDown.Value = currentParam;
                else
                    app.ParamDropDown.Value = common{1};
                end
                app.StatusLabel.Text = sprintf('Ready: %d common parameters found', numel(common));
            else
                app.ParamDropDown.Items = {'No common parameters'};
                app.ParamDropDown.Value = 'No common parameters';
                app.StatusLabel.Text = 'Warning: No matching parameters between ABL and CDI';
            end
        end

        % --- DATA INGESTION: PARSE ABL ---
        function [tbl, patientIDs] = parseABL(~, fullpath)
            fid = fopen(fullpath, 'r');
            if fid == -1
                error('ABL file: cannot open file.');
            end
            firstLine = fgetl(fid);
            fclose(fid);

            if length(strfind(firstLine, ';')) > length(strfind(firstLine, ','))
                delimiter = ';';
            else
                delimiter = ',';
            end

            numCols = length(strfind(firstLine, delimiter)) + 1;
            formatSpec = repmat('%q', 1, numCols);
            fid = fopen(fullpath, 'r');
            dataCell = textscan(fid, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);
            fclose(fid);

            numRows = length(dataCell{1});
            raw = cell(numRows, numCols);
            for col = 1 : numCols
                for row = 1 : numRows
                    if col <= length(dataCell) && row <= length(dataCell{col})
                        raw{row, col} = dataCell{col}{row};
                    else
                        raw{row, col} = '';
                    end
                end
            end

            headers = raw(1, :);
            timeStrings = strtrim(raw(2:end, 1));
            timeVec     = datetime(timeStrings, 'InputFormat', 'd.M.yyyy HH:mm');
            
            patientIDCol = find(strcmpi(strtrim(headers), 'Patient Id') | strcmpi(strtrim(headers), 'PatientId') | strcmpi(strtrim(headers), 'Patient_Id'));
            if isempty(patientIDCol)
                for k = 1:numel(headers)
                    if contains(lower(strtrim(headers{k})), 'patient') && contains(lower(strtrim(headers{k})), 'id')
                        patientIDCol = k;
                        break;
                    end
                end
            end
            
            if ~isempty(patientIDCol)
                patientIDValues = strtrim(raw(2:end, patientIDCol(1)));
                patientIDValues = strrep(patientIDValues, '"', '');
            else
                patientIDValues = repmat({''}, numel(timeStrings), 1);
            end

            keepIdx  = [];
            keepName = {};
            keepData = {}; 
            for k = 2 : numCols
                hdr = strtrim(headers{k});
                if strcmpi(hdr, 'Error code') || strcmpi(hdr, 'Error_code') || isempty(hdr)
                    continue;
                end
                if contains(lower(hdr), 'patient') || contains(lower(hdr), 'last name') || ...
                   contains(lower(hdr), 'first name') || contains(lower(hdr), 'sample #') || ...
                   contains(lower(hdr), 'status') || contains(lower(hdr), 'sample type') || ...
                   contains(lower(hdr), 'sex') || contains(lower(hdr), 'birthdate') || ...
                   contains(lower(hdr), 'note') || contains(lower(hdr), 'department') || ...
                   contains(lower(hdr), 'accession') || contains(lower(hdr), 'site') || ...
                   contains(lower(hdr), 'draw time') || contains(lower(hdr), 'physician') || ...
                   contains(lower(hdr), 'operator') || contains(lower(hdr), 'approval') || ...
                   contains(lower(hdr), 'report layout') || contains(lower(hdr), 'measuring mode') || ...
                   contains(lower(hdr), 'errors detected')
                    continue;
                end
                col     = strtrim(raw(2:end, k));
                col     = strrep(col, ',', '.');
                numVals = str2double(col);
                if any(~isnan(numVals))
                    keepIdx(end+1)  = k;               
                    cleanName = regexprep(hdr, '\s*\(.*?\)\s*', '');
                    cleanName = strrep(cleanName, '"', '');
                    cleanName = strtrim(cleanName);
                    keepName{end+1} = cleanName;       
                    keepData{end+1} = numVals;         
                end
            end

            if isempty(keepIdx)
                error('ABL file: no usable numeric columns found.');
            end

            seen = containers.Map('KeyType','char','ValueType','double');
            seenIdx = containers.Map('KeyType','char','ValueType','double');
            toRemove = false(1, numel(keepName));
            
            for n = 1 : numel(keepName)
                nm = keepName{n};
                currentData = keepData{n};
                currentScore = sum(~isnan(currentData) & currentData ~= 0);
                
                if seen.isKey(nm)
                    prevScore = seen(nm);
                    prevIdx = seenIdx(nm);
                    
                    if currentScore > prevScore
                        toRemove(prevIdx) = true;
                        seen(nm) = currentScore;
                        seenIdx(nm) = n;
                    else
                        toRemove(n) = true;
                    end
                else
                    seen(nm) = currentScore;
                    seenIdx(nm) = n;
                end
            end
            
            keepIdx = keepIdx(~toRemove);
            keepName = keepName(~toRemove);
            keepData = keepData(~toRemove);

            tbl = table('Size',          [numel(timeVec), 2+numel(keepIdx)], ...
                        'VariableTypes', [{'datetime'}, {'cell'}, repmat({'double'}, 1, numel(keepIdx))], ...
                        'VariableNames', [{'Time'}, {'PatientID'}, keepName]);

            tbl.Time = timeVec;
            tbl.PatientID = patientIDValues;
            for k = 1 : numel(keepIdx)
                tbl{:, k+2} = keepData{k};
            end

            tbl = tbl(~isnat(tbl.Time), :);
            tbl = sortrows(tbl, 'Time');
            
            uniqueIDs = unique(tbl.PatientID);
            uniqueIDs = uniqueIDs(~cellfun(@isempty, uniqueIDs));
            numericIDs = str2double(uniqueIDs);
            if all(~isnan(numericIDs))
                [~, sortIdx] = sort(numericIDs);
                patientIDs = uniqueIDs(sortIdx);
            else
                patientIDs = sort(uniqueIDs);
            end
        end

        % --- DATA INGESTION: PARSE CDI ---
        function tbl = parseCDI(~, fullpath)
            allColNames = {'Time','pH','pCO2','pO2','TEMP','HCO3','BE','sO2','K+',...
                        'VO2','Q','BSA','pH_v','pCO2_v','pO2_v','TEMP_v','SO2_v','HCT','tHb'};

            fid = fopen(fullpath, 'r', 'native', 'UTF-8');
            if fid == -1
                error('CDI file: cannot open file.');
            end
            rawText = fread(fid, Inf, 'char=>char')';
            fclose(fid);

            if length(rawText) >= 3 && ...
               rawText(1) == char(239) && rawText(2) == char(187) && rawText(3) == char(191)
                rawText = rawText(4:end);
            end

            rawText   = strrep(rawText, char(13), '');
            lines     = strsplit(rawText, char(10));
            lines     = strtrim(lines);
            lines(cellfun(@isempty, lines)) = [];   

            nLines = numel(lines);
            
            delimiter = ','; 
            for i = 1:min(20, nLines)
                ln = lines{i};
                if ~isempty(ln) && ln(1) == '['
                    numTabs = sum(ln == char(9));
                    numCommas = sum(ln == ',');
                    if numTabs > numCommas
                        delimiter = char(9); 
                    end
                    break;
                end
            end

            timeOut  = NaT(nLines, 1);
            dataOut  = nan(nLines, numel(allColNames)-1);   
            validRow = false(nLines, 1);

            for i = 1 : nLines
                ln = lines{i};

                if isempty(ln) || ln(1) ~= '['
                    continue;
                end
                bracketEnd = strfind(ln, ']');
                if isempty(bracketEnd)
                    continue;
                end
                bracketEnd = bracketEnd(1);

                logDateStr = ln(2 : bracketEnd-1);            
                remainder  = strtrim(ln(bracketEnd+1 : end)); 

                if isempty(remainder) || isletter(remainder(1))
                    continue;
                end

                fields = strsplit(remainder, delimiter);
                fields = strtrim(fields);

                devTime = fields{1};
                if sum(devTime == ':') ~= 2 || length(devTime) > 8
                    continue;
                end

                try
                    logDate   = datetime(logDateStr, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
                    fullStamp = datetime([datestr(logDate,'yyyy-mm-dd') ' ' devTime], ...
                                         'InputFormat', 'yyyy-MM-dd HH:mm:ss');
                catch
                    continue;
                end

                if numel(fields) < 2
                    continue;
                end

                nFields  = min(numel(fields), numel(allColNames));
                rowVals  = nan(1, numel(allColNames)-1);
                anyValid = false;
                for f = 2 : nFields
                    val = strtrim(fields{f});
                    val = strrep(val, '%', '');        
                    val = strrep(val, '---', '');      
                    val = strrep(val, '--', '');       
                    val = strrep(val, '-.-', '');      
                    num = str2double(val);
                    if ~isnan(num)
                        rowVals(f-1) = num;
                        anyValid     = true;
                    end
                end

                if ~anyValid
                    continue;
                end

                timeOut(i)   = fullStamp;
                dataOut(i,:) = rowVals;
                validRow(i)  = true;
            end

            timeOut  = timeOut(validRow);
            dataOut  = dataOut(validRow, :);

            tbl = table(timeOut, 'VariableNames', {'Time'});
            
            for c = 2 : numel(allColNames)
                colData = dataOut(:, c-1);
                if any(~isnan(colData))  
                    tbl.(allColNames{c}) = colData;
                end
            end

            tbl = sortrows(tbl, 'Time');
        end

        function TimeShiftSpinnerValueChanged(app, ~)
            if ~isempty(app.ABL_Table) && ~isempty(app.CDI_Table)
                app.StatusLabel.Text = 'Recalculating shift...'; 
                drawnow;
                AnalyzeButtonPushed(app, []);
            end
        end

        function FitWindowCheckBoxChanged(app, ~)
            if ~isempty(app.ABL_Table) && ~isempty(app.CDI_Table)
                drawFitWindowShade(app);
                AnalyzeButtonPushed(app, []);
            end
        end

        function FitWindowApplyButtonPushed(app, ~)
            if isempty(app.ABL_Table) || isempty(app.CDI_Table)
                return;
            end

            startStr = strtrim(app.FitWindowStartEdit.Value);
            endStr   = strtrim(app.FitWindowEndEdit.Value);

            if isempty(startStr) && isempty(endStr)
                if isfield(app.Stats, 'overlapStart') && ~isnat(app.Stats.overlapStart)
                    app.FitWindowStartEdit.Value = datestr(app.Stats.overlapStart, 'dd.mm.yyyy HH:MM');
                    app.FitWindowEndEdit.Value   = datestr(app.Stats.overlapEnd,   'dd.mm.yyyy HH:MM');
                end
            else
                ok = true;
                try
                    datetime(startStr, 'InputFormat', 'dd.MM.yyyy HH:mm');
                catch
                    ok = false;
                    uialert(app.UIFigure, ...
                        sprintf('Start time format invalid: "%s"\nUse: dd.MM.yyyy HH:mm  (e.g. 30.01.2026 11:00)', startStr), ...
                        'Invalid Window Start');
                end
                if ok
                    try
                        datetime(endStr, 'InputFormat', 'dd.MM.yyyy HH:mm');
                    catch
                        ok = false;
                        uialert(app.UIFigure, ...
                            sprintf('End time format invalid: "%s"\nUse: dd.MM.yyyy HH:mm  (e.g. 30.01.2026 16:00)', endStr), ...
                            'Invalid Window End');
                    end
                end
                if ~ok, return; end
            end

            app.FitWindowCheckBox.Value = true;

            drawFitWindowShade(app);
            AnalyzeButtonPushed(app, []);
        end

        function TimeToleranceSpinnerValueChanged(app, ~)
            if ~isempty(app.ABL_Table) && ~isempty(app.CDI_Table)
                app.StatusLabel.Text = 'Recalculating tolerance...'; 
                drawnow;
                AnalyzeButtonPushed(app, []);
            end
        end

        function FitWindowAutoButtonPushed(app, ~)
            if isempty(app.ABL_Table) || isempty(app.CDI_Table)
                uialert(app.UIFigure, 'Please load both ABL and CDI data first.', 'No Data');
                return;
            end

            param = app.ParamDropDown.Value;

            shiftMins = app.TimeShiftSpinner.Value;
            cdi_times = app.CDI_Table.Time + minutes(shiftMins);
            if ismember(param, app.CDI_Table.Properties.VariableNames)
                cdi_vals = app.CDI_Table.(param);
            else
                cdi_vals = ones(height(app.CDI_Table), 1);
            end
            validCDI  = ~isnan(cdi_vals);
            cdi_start = min(cdi_times(validCDI));
            cdi_end   = max(cdi_times(validCDI));

            abl_times = app.ABL_Table.Time;
            if ismember(param, app.ABL_Table.Properties.VariableNames)
                abl_vals = app.ABL_Table.(param);
            else
                abl_vals = ones(height(app.ABL_Table), 1);
            end
            validABL = ~isnan(abl_vals);
            abl_times = abl_times(validABL);

            abl_in_cdi = abl_times(abl_times >= cdi_start & abl_times <= cdi_end);

            if isempty(abl_in_cdi)
                uialert(app.UIFigure, ...
                    'No ABL draws found within the CDI recording window. Check time shift.', ...
                    'No Overlap');
                return;
            end

            autoStart = min(abl_in_cdi);
            autoEnd   = max(abl_in_cdi);

            app.FitWindowStartEdit.Value = datestr(autoStart, 'dd.mm.yyyy HH:MM');
            app.FitWindowEndEdit.Value   = datestr(autoEnd,   'dd.mm.yyyy HH:MM');

            app.FitWindowCheckBox.Value = true;
            drawFitWindowShade(app);
            AnalyzeButtonPushed(app, []);

            app.StatusLabel.Text = sprintf('Auto window: %s → %s (%d ABL draws)', ...
                datestr(autoStart,'HH:MM'), datestr(autoEnd,'HH:MM'), numel(abl_in_cdi));
        end

        % --- AUTO TIME SHIFT DETECTION ---
        function AutoShiftButtonPushed(app, ~)
            if isempty(app.ABL_Table) || isempty(app.CDI_Table)
                uialert(app.UIFigure, 'Please load both ABL and CDI data files first.', 'Missing Data');
                return;
            end
            param = app.ParamDropDown.Value;
            if isempty(param) || strcmpi(param, 'Load files first...')
                uialert(app.UIFigure, 'Please select a valid parameter first.', 'No Parameter');
                return;
            end
            
            ablForAlign = app.ABL_Table;
            if ismember('PatientID', ablForAlign.Properties.VariableNames)
                ablForAlign = removevars(ablForAlign, 'PatientID');
            end
            abl_tt = table2timetable(ablForAlign);
            base_cdi_tt = table2timetable(app.CDI_Table);
            
            vA = abl_tt.Properties.VariableNames;
            vC = base_cdi_tt.Properties.VariableNames;
            aName = ''; cName = '';
            for i=1:length(vA), if strcmpi(vA{i}, param), aName=vA{i}; break; end; end
            for i=1:length(vC), if strcmpi(vC{i}, param), cName=vC{i}; break; end; end
            
            if isempty(aName) || isempty(cName)
                uialert(app.UIFigure, 'Column missing in one of the files.', 'Data Error');
                app.StatusLabel.Text = 'Auto-shift failed.';
                return;
            end
            
            tA = abl_tt.Time;
            yA = abl_tt.(aName);
            tC = base_cdi_tt.Time;
            yC = base_cdi_tt.(cName);
            
            validA = ~isnan(yA);
            tA = tA(validA); yA = yA(validA);
            validC = ~isnan(yC);
            tC = tC(validC); yC = yC(validC);

            useFitWin = false;
            fitWinStart = NaT; fitWinEnd = NaT;
            if app.FitWindowCheckBox.Value && ...
                    ~isempty(strtrim(app.FitWindowStartEdit.Value))
                try
                    fitWinStart = datetime(app.FitWindowStartEdit.Value, ...
                        'InputFormat', 'dd.MM.yyyy HH:mm');
                    fitWinEnd   = datetime(app.FitWindowEndEdit.Value, ...
                        'InputFormat', 'dd.MM.yyyy HH:mm');
                    useFitWin = true;
                catch
                end
            end

            if useFitWin
                winMaskA = tA >= fitWinStart & tA <= fitWinEnd;
                tA = tA(winMaskA); yA = yA(winMaskA);
                if numel(tA) < 4
                    uialert(app.UIFigure, ...
                        'Fitting window contains fewer than 4 ABL draws. Widen the window or disable it for Auto-Shift.', ...
                        'Window Too Narrow');
                    app.StatusLabel.Text = 'Auto-shift failed - window too narrow.';
                    return;
                end
                app.StatusLabel.Text = sprintf('Running shift search inside fitting window (%d ABL draws)...', numel(tA)); drawnow;
            else
                app.StatusLabel.Text = 'Running simulated shifts (-120 to +120 min)...'; drawnow;
            end
            
            tolDays = minutes(app.TimeToleranceSpinner.Value);
            shifts = -120:1:120;
            r_vals = nan(size(shifts));
            n_vals = zeros(size(shifts));
            
            for i = 1:length(shifts)
                shifted_tC = tC + minutes(shifts(i));
                tC_search  = shifted_tC;
                yC_search  = yC;

                if useFitWin
                    winMaskC  = tC_search >= fitWinStart & tC_search <= fitWinEnd;
                    tC_search = tC_search(winMaskC);
                    yC_search = yC_search(winMaskC);
                end
                if isempty(tC_search), continue; end

                mA = []; mC = [];
                for j = 1:length(tA)
                    [md, idx] = min(abs(tC_search - tA(j)));
                    if md <= tolDays
                        mA(end+1) = yA(j); 
                        mC(end+1) = yC_search(idx); 
                    end
                end
                
                n_vals(i) = length(mA);
                if n_vals(i) >= 4 && std(mA) > 0 && std(mC) > 0
                    R = corrcoef(mA, mC);
                    r_vals(i) = R(1,2);
                end
            end
            
            max_possible_pairs = max(n_vals);
            valid_idx = n_vals >= max_possible_pairs * 0.5 & n_vals >= 4;
            
            valid_shifts = shifts(valid_idx);
            valid_r = r_vals(valid_idx);
            
            if isempty(valid_r) || all(isnan(valid_r))
                uialert(app.UIFigure, 'Could not find a valid pattern match. The data might be flatlined or disconnected.', 'Auto-Shift Failed');
                app.StatusLabel.Text = 'Auto-shift failed.';
                return;
            end
            
            [max_r, best_idx] = max(valid_r);
            best_shift = valid_shifts(best_idx);
            
            app.TimeShiftSpinner.Value = best_shift;
            
            if useFitWin
                app.StatusLabel.Text = sprintf('Best correlation shift (window): %d min (r = %.3f)', best_shift, max_r);
            else
                app.StatusLabel.Text = sprintf('Best correlation shift: %d min (r = %.3f)', best_shift, max_r);
            end
            AnalyzeButtonPushed(app, []); 
        end
        
        % --- AUTO-TUNE HYBRID & WEIGHTED DEMING PARAMETERS ---
        function AutoTuneButtonPushed(app, ~)
            if isempty(app.Aligned_Data) || ~isfield(app.Stats, 'TimeShift')
                uialert(app.UIFigure, 'Please run Analysis first to extract paired data.', 'No Data');
                return;
            end
            
            param = app.CurrentParam;
            ablCol = ''; cdiCol = '';
            varNames = app.Aligned_Data.Properties.VariableNames;
            for i = 1:length(varNames)
                v = varNames{i};
                if contains(v, param) && contains(v, 'abl'), ablCol = v; end
                if contains(v, param) && contains(v, 'cdi'), cdiCol = v; end
            end
            if isempty(ablCol) || isempty(cdiCol)
                uialert(app.UIFigure, 'Could not find valid columns for tuning.', 'Error');
                return; 
            end
            
            xABL = app.Aligned_Data.(ablCol);
            yCDI = app.Aligned_Data.(cdiCol);
            validIdx = ~isnan(xABL) & ~isnan(yCDI);
            xABL = xABL(validIdx);
            validTimes = app.Aligned_Data.Time(validIdx);
            
            winMask = getFitWindowMask(app, validTimes);
            xABL = xABL(winMask);
            validTimes = validTimes(winMask);
            
            if numel(xABL) < 3
                uialert(app.UIFigure, 'Not enough data points in the window for accurate tuning.', 'Insufficient Data');
                return;
            end
            
            app.StatusLabel.Text = 'Auto-tuning Hybrid & Weighted Deming Parameters...'; drawnow;
            
            shiftMins = app.Stats.TimeShift;
            fullCDITime = app.CDI_Table.Time + minutes(shiftMins);
            fullCDIVals = app.CDI_Table.(param);
            
            validCDIMask = ~isnan(fullCDIVals);
            searchTimes = fullCDITime(validCDIMask);
            matchIndices = zeros(size(xABL));
            for k = 1:numel(xABL)
                [~, bestIdx] = min(abs(searchTimes - validTimes(k)));
                matchIndices(k) = bestIdx;
            end
            
            raw_cdi_paired = fullCDIVals(validCDIMask);
            raw_cdi_paired = raw_cdi_paired(matchIndices);
            
            globalCleanMask = robustCleanMask(app, xABL, raw_cdi_paired);
            xClean_global = xABL(globalCleanMask);
            
            best_rmse = inf;
            best_tau_r = 0.0;
            best_tau_f = 0.0;
            best_lam = 1.0;
            best_w1  = app.SmoothW1Spinner.Value;
            
            tau_cands = 0:1.0:8;
            lam_cands = [0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0];
            w1_cands  = [4, 8, 12, 16, 24, 32];
            
            winMaskCDI = searchTimes >= min(validTimes) & searchTimes <= max(validTimes);
            raw_valid_full = fullCDIVals(validCDIMask);
            raw_valid_win = raw_valid_full(winMaskCDI);
            
            if numel(raw_valid_win) > 2
                raw_roughness = std(diff(raw_valid_win(~isnan(raw_valid_win))), 'omitnan');
                if raw_roughness == 0 || isnan(raw_roughness), raw_roughness = 1; end
            else
                raw_roughness = 1;
            end
            
            for w1 = w1_cands
                for tr = tau_cands
                    tf = tr; % Strict alignment with LOO-CV symmetric rise=fall rule
                    cdi_fast_full = computeAsymmetricFastCDI(app, fullCDIVals, fullCDITime, w1, tr, tf);
                    searchVals = cdi_fast_full(validCDIMask);
                    cdi_fast_paired = searchVals(matchIndices);
                    
                    sv_win = searchVals(winMaskCDI);
                    rv = sv_win(~isnan(sv_win));
                    if numel(rv) > 2
                        base_corr_rough = std(diff(rv), 'omitnan');
                    else
                        base_corr_rough = raw_roughness;
                    end
                    
                    for l = lam_cands
                        [slope, intercept] = fitWeightedDeming(app, xClean_global, cdi_fast_paired(globalCleanMask), l);
                        yCorr = (cdi_fast_paired - intercept) / slope;
                        rmse = sqrt(mean((yCorr(globalCleanMask) - xClean_global).^2, 'omitnan'));
                        
                        corr_roughness = base_corr_rough / abs(slope);
                        if isnan(corr_roughness), corr_roughness = raw_roughness; end
                        roughness_ratio = corr_roughness / raw_roughness;
                        
                        smoothness_penalty = max(0, roughness_ratio - 1.0) * 0.25 * rmse;
                        penalised_rmse = rmse + smoothness_penalty;
                        
                        if penalised_rmse < best_rmse
                            best_rmse = penalised_rmse;
                            best_tau_r = tr;
                            best_tau_f = tf;
                            best_lam = l;
                            best_w1  = w1;
                        end
                    end
                end
            end
            
            app.TauSpinner.Value = best_tau_r;
            app.TauFallSpinner.Value = best_tau_f;
            app.SmoothW1Spinner.Value = best_w1;
            app.DemingLambdaEditField.Value = best_lam;
            app.CorrectionMethodDropDown.Value = 'Hybrid (Time-Series + Deming)';
            app.StatusLabel.Text = sprintf('Auto-Tune: τ_r=%.1f τ_f=%.1f | W1=%d | λ=%.2f', best_tau_r, best_tau_f, best_w1, best_lam);
            
            ApplyCorrectionButtonPushed(app, []);
        end

        function CorrViewSwitchChanged(app, ~)
            if isempty(app.Stats) || ~isfield(app.Stats, 'xABL')
                return;
            end

            param = app.CurrentParam;
            xABL  = app.Stats.xABL;

            cla(app.CorrelationAxes);
            hold(app.CorrelationAxes, 'on');

            if app.CorrViewSwitch.Value && isfield(app.CorrectionModel, 'yCorrected') && ~isempty(app.CorrectionModel.yCorrected)
                yVals = app.CorrectionModel.yCorrected;
                if numel(yVals) ~= numel(xABL)
                    n_min = min(numel(yVals), numel(xABL));
                    yVals = yVals(1:n_min); xABL = xABL(1:n_min);
                end
                p_val = app.CorrectionModel.pCorr;
                r_val = app.CorrectionModel.r_new;
                sc = scatter(app.CorrelationAxes, xABL, yVals, 50, 'filled', ...
                    'MarkerFaceColor', [0.2 0.7 0.2], 'DisplayName', sprintf('Corrected (r=%.3f)', r_val));
                sc.DataTipTemplate.DataTipRows(1).Label = 'ABL';
                sc.DataTipTemplate.DataTipRows(2).Label = 'CDI (Corr)';
                yLabelStr = ['Corrected CDI ' param];
                titleStr  = sprintf('After Correction (r=%.3f)', r_val);
            else
                yVals = app.Stats.yCDI;
                % Mask baseline to match the currently applied MAD filter
                if isfield(app.CorrectionModel, 'yCorrected') && ~isempty(app.CorrectionModel.yCorrected)
                    yCorrTmp = app.CorrectionModel.yCorrected;
                    if numel(yCorrTmp) == numel(xABL)
                        validPairs = ~isnan(yCorrTmp) & ~isnan(xABL);
                        yVals = yVals(validPairs);
                        xABL = xABL(validPairs);
                    end
                end
                if numel(xABL) >= 2 && std(xABL, 'omitnan') > 0
                    p_val = polyfit(xABL, yVals, 1);
                    R_tmp = corrcoef(xABL, yVals);
                    r_val = R_tmp(1,2);
                else
                    p_val = app.Stats.p; r_val = app.Stats.r;
                end
                sc = scatter(app.CorrelationAxes, xABL, yVals, 50, 'filled', ...
                    'MarkerFaceColor', [0 0.4470 0.7410], 'DisplayName', sprintf('Paired samples (r=%.3f)', r_val));
                sc.DataTipTemplate.DataTipRows(1).Label = 'ABL';
                sc.DataTipTemplate.DataTipRows(2).Label = 'CDI';
                yLabelStr = ['CDI ' param];
                titleStr  = sprintf('Correlation (r=%.3f)', r_val);
            end

            validX = xABL(~isnan(xABL) & ~isinf(xABL));
            validY = yVals(~isnan(yVals) & ~isinf(yVals));
            if isempty(validX) || isempty(validY)
                minV = 0; maxV = 1;
            else
                minV = min([validX; validY]); maxV = max([validX; validY]);
            end
            span = maxV - minV;
            if span == 0 || isnan(span), span = 0.1; end
            minV = minV - span * 0.1; maxV = maxV + span * 0.1;
            if minV >= maxV || isnan(minV) || isnan(maxV), minV = 0; maxV = 1; end

            plot(app.CorrelationAxes, [minV maxV], [minV maxV], 'k--', 'DisplayName', 'Identity');
            if ~any(isnan(p_val))
                xFit = linspace(minV, maxV, 100);
                plot(app.CorrelationAxes, xFit, polyval(p_val, xFit), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Regression');
            end

            xlabel(app.CorrelationAxes, ['ABL ' param]);
            ylabel(app.CorrelationAxes, yLabelStr);
            title(app.CorrelationAxes, titleStr);
            legend(app.CorrelationAxes, 'Location', 'best');
            axis(app.CorrelationAxes, 'square');
            xlim(app.CorrelationAxes, [minV, maxV]);
            ylim(app.CorrelationAxes, [minV, maxV]);
            grid(app.CorrelationAxes, 'on');
            hold(app.CorrelationAxes, 'off');
        end

        function BAViewSwitchChanged(app, ~)
            if isempty(app.Stats) || ~isfield(app.Stats, 'xABL')
                return;
            end

            xABL = app.Stats.xABL;
            cla(app.BlandAltmanAxes);
            hold(app.BlandAltmanAxes, 'on');

            if app.BAViewSwitch.Value && isfield(app.CorrectionModel, 'yCorrected') && ~isempty(app.CorrectionModel.yCorrected)
                yVals = app.CorrectionModel.yCorrected;
                if numel(yVals) ~= numel(xABL)
                    n_min = min(numel(yVals), numel(xABL));
                    yVals = yVals(1:n_min); xABL = xABL(1:n_min);
                end
                diffVals = yVals - xABL;
                bias = mean(diffVals, 'omitnan'); sd = std(diffVals, 'omitnan');
                methodLabel = app.CorrectionMethodDropDown.Value;
                if isfield(app.CorrectionModel, 'autoSelected') && app.CorrectionModel.autoSelected
                    methodLabel = app.CorrectionModel.autoWinner;
                end
                titleStr  = sprintf('After %s', methodLabel);
                markerColor = [0.2 0.7 0.2];
                legendName  = sprintf('Corrected  Bias=%.4f  SD=%.4f', bias, sd);
                yLabelStr   = 'Diff (CDI corr - ABL)';
                xLabelStr   = 'Mean (ABL, CDI corrected)';
                lineColor   = [0.1 0.6 0.1];
            else
                yVals = app.Stats.yCDI;
                % Mask baseline to match the currently applied MAD filter
                if isfield(app.CorrectionModel, 'yCorrected') && ~isempty(app.CorrectionModel.yCorrected)
                    yCorrTmp = app.CorrectionModel.yCorrected;
                    if numel(yCorrTmp) == numel(xABL)
                        validPairs = ~isnan(yCorrTmp) & ~isnan(xABL);
                        yVals = yVals(validPairs);
                        xABL = xABL(validPairs);
                    end
                end
                diffVals = yVals - xABL;
                bias = mean(diffVals, 'omitnan'); sd = std(diffVals, 'omitnan');
                titleStr    = 'Bland-Altman';
                markerColor = [0 0.4470 0.7410];
                legendName  = sprintf('Data  Bias=%.4f  SD=%.4f', bias, sd);
                yLabelStr   = 'Diff (CDI - ABL)';
                xLabelStr   = 'Mean (ABL, CDI)';
                lineColor   = [0 0.4470 0.7410];
            end

            avgVals = (xABL + yVals) / 2;
            validAvg = avgVals(~isnan(avgVals) & ~isinf(avgVals));
            if isempty(validAvg)
                xlims = [0 1];
            else
                span = max(validAvg) - min(validAvg);
                if span == 0 || isnan(span), span = 0.1; end
                xlims = [min(validAvg)-span*0.2, max(validAvg)+span*0.2];
                if xlims(1) >= xlims(2) || isnan(xlims(1)), xlims = [0 1]; end
            end

            loa_up = bias + 1.96*sd; loa_lo = bias - 1.96*sd;

            sc = scatter(app.BlandAltmanAxes, avgVals, diffVals, 50, 'filled', ...
                'MarkerFaceColor', markerColor, 'DisplayName', legendName);
            sc.DataTipTemplate.DataTipRows(1).Label = 'Mean';
            sc.DataTipTemplate.DataTipRows(2).Label = 'Diff';

            if ~isnan(bias)
                plot(app.BlandAltmanAxes, xlims, [bias bias], '-', 'Color', lineColor, ...
                    'LineWidth', 1.5, 'DisplayName', sprintf('Bias = %.4f', bias));
                plot(app.BlandAltmanAxes, xlims, [loa_up loa_up], '--', 'Color', lineColor, ...
                    'LineWidth', 1, 'DisplayName', sprintf('+1.96SD = %.4f', loa_up));
                plot(app.BlandAltmanAxes, xlims, [loa_lo loa_lo], '--', 'Color', lineColor, ...
                    'LineWidth', 1, 'DisplayName', sprintf('-1.96SD = %.4f', loa_lo));
                plot(app.BlandAltmanAxes, xlims, [0 0], 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
            end

            xlabel(app.BlandAltmanAxes, xLabelStr);
            ylabel(app.BlandAltmanAxes, yLabelStr);
            title(app.BlandAltmanAxes, titleStr);
            xlim(app.BlandAltmanAxes, xlims);
            legend(app.BlandAltmanAxes, 'Location', 'best');
            grid(app.BlandAltmanAxes, 'on');
            hold(app.BlandAltmanAxes, 'off');
        end

        % --- CORE ANALYSIS AND SYNCHRONIZATION ---
        function AnalyzeButtonPushed(app, ~)
            if isempty(app.ABL_Table) || isempty(app.CDI_Table)
                uialert(app.UIFigure, 'Please load both ABL and CDI data files first.', 'Missing Data');
                return;
            end
            param = app.ParamDropDown.Value;
            if isempty(param) || strcmpi(param, 'Load files first...')
                uialert(app.UIFigure, 'Please select a valid parameter first.', 'No Parameter');
                return;
            end
            app.CurrentParam = param;

            ablCols = setdiff(app.ABL_Table.Properties.VariableNames, {'PatientID'});
            if ~ismember(param, ablCols)
                uialert(app.UIFigure, ['Parameter "' param '" not found in ABL data.'], 'Missing Parameter');
                return;
            end
            if ~ismember(param, app.CDI_Table.Properties.VariableNames)
                uialert(app.UIFigure, ['Parameter "' param '" not found in CDI data.'], 'Missing Parameter');
                return;
            end

            app.StatusLabel.Text = 'Analyzing...'; drawnow;

            tolMins = app.TimeToleranceSpinner.Value;

            ablForAlign = app.ABL_Table;
            if ismember('PatientID', ablForAlign.Properties.VariableNames)
                ablForAlign = removevars(ablForAlign, 'PatientID');
            end

            abl_tt = table2timetable(ablForAlign);
            
            shiftMins = app.TimeShiftSpinner.Value;
            shifted_CDI = app.CDI_Table;
            shifted_CDI.Time = shifted_CDI.Time + minutes(shiftMins);
            cdi_tt = table2timetable(shifted_CDI);

            aligned = synchronize(abl_tt, cdi_tt, 'first', 'nearest');

            cdi_times = cdi_tt.Time;
            n = height(aligned);
            timeDiff = nan(n,1);
            for k = 1:n
                abl_t = aligned.Time(k);
                [minD, ~] = min(abs(cdi_times - abl_t));
                timeDiff(k) = minutes(minD);
            end

            validMask = timeDiff <= tolMins;
            aligned = aligned(validMask, :);

            if height(aligned) == 0
                uialert(app.UIFigure, 'No paired samples within tolerance. Try increasing tolerance.', 'No Pairs');
                app.StatusLabel.Text = 'No pairs found.';
                return;
            end

            varNames  = aligned.Properties.VariableNames;
            ablCol = '';
            cdiCol = '';
            for i=1:length(varNames)
                v = varNames{i};
                if endsWith(v, '_abl_tt') && startsWith(v, param)
                    ablCol = v;
                end
                if endsWith(v, '_cdi_tt') && startsWith(v, param)
                    cdiCol = v;
                end
            end
            if isempty(ablCol)
                for i=1:length(varNames)
                    v = varNames{i};
                    if contains(v, param) && contains(v, 'abl')
                        ablCol = v; break;
                    end
                end
            end
            if isempty(cdiCol)
                for i=1:length(varNames)
                    v = varNames{i};
                    if contains(v, param) && contains(v, 'cdi')
                        cdiCol = v; break;
                    end
                end
            end
            if isempty(ablCol) || isempty(cdiCol)
                uialert(app.UIFigure, ['Could not find matched columns for "' param '".'], 'Column Error');
                app.StatusLabel.Text = 'Column match failed.';
                return;
            end

            xABL = aligned.(ablCol);
            yCDI = aligned.(cdiCol);
            validIdx = ~isnan(xABL) & ~isnan(yCDI);
            xABL = xABL(validIdx);
            yCDI = yCDI(validIdx);
            timeVals = aligned.Time(validIdx);

            if isempty(xABL)
                uialert(app.UIFigure, 'No valid data pairs after filtering NaN values.', 'No Valid Data');
                app.StatusLabel.Text = 'No valid pairs.';
                return;
            end

            app.Aligned_Data = aligned(validIdx, :);

            statsMask = getFitWindowMask(app, timeVals);
            xABL_stat = xABL(statsMask);
            yCDI_stat = yCDI(statsMask);

            if numel(xABL_stat) < 2
                xABL_stat = xABL;
                yCDI_stat = yCDI;
                if app.FitWindowCheckBox.Value
                    app.StatusLabel.Text = 'Warning: fitting window has <2 pairs - using all pairs for statistics.';
                end
            end

            diff_vals = yCDI_stat - xABL_stat;
            bias = mean(diff_vals, 'omitnan');
            sd   = std(diff_vals, 'omitnan');
            loa_up = bias + 1.96*sd;
            loa_lo = bias - 1.96*sd;
            if std(xABL_stat,'omitnan') > 0 && std(yCDI_stat,'omitnan') > 0
                R = corrcoef(xABL_stat, yCDI_stat);
                r = R(1,2);
            else
                r = NaN;
            end

            p = polyfit(xABL_stat, yCDI_stat, 1);
            slope = p(1);
            intercept = p(2);

            raw_abl_t_pre = abl_tt.Time;
            raw_abl_y_pre = abl_tt.(param);
            full_cdi_t_pre = shifted_CDI.Time;
            full_cdi_y_pre = shifted_CDI.(param);
            valid_abl_pre = ~isnan(raw_abl_y_pre);
            valid_cdi_pre = ~isnan(full_cdi_y_pre);

            cdi_start = min(full_cdi_t_pre(valid_cdi_pre));
            cdi_end   = max(full_cdi_t_pre(valid_cdi_pre));
            abl_start = min(raw_abl_t_pre(valid_abl_pre));
            abl_end   = max(raw_abl_t_pre(valid_abl_pre));

            overlap_start = max(cdi_start, abl_start);
            overlap_end   = min(cdi_end,   abl_end);
            hasOverlap    = overlap_start < overlap_end;

            app.Stats = struct('N', numel(xABL_stat), 'Ntotal', numel(xABL), ...
                               'bias', bias, 'sd', sd, ...
                               'loa_up', loa_up, 'loa_lo', loa_lo, 'r', r, ...
                               'slope', slope, 'intercept', intercept, ...
                               'param', param, ...
                               'xABL', xABL_stat, 'yCDI', yCDI_stat, 'p', p, ...
                               'TimeShift', shiftMins, ...
                               'overlapStart', overlap_start, ...
                               'overlapEnd', overlap_end);

            if app.FitWindowCheckBox.Value && numel(xABL_stat) < numel(xABL)
                app.NPairsLabel.Text = sprintf('N Pairs: %d (window) / %d total', numel(xABL_stat), numel(xABL));
            else
                app.NPairsLabel.Text = sprintf('N Pairs: %d', app.Stats.N);
            end
            app.BiasLabel.Text   = sprintf('Bias: %.3f', app.Stats.bias);
            app.SDLabel.Text     = sprintf('SD: %.3f', app.Stats.sd);
            app.LOALabel.Text    = sprintf('95%% LoA: [%.3f, %.3f]', app.Stats.loa_lo, app.Stats.loa_up);
            app.CorrelationLabel.Text = sprintf('r = %.4f', app.Stats.r);
            app.RegressionLabel.Text  = sprintf('y = %.3fx + %.3f', app.Stats.slope, app.Stats.intercept);
            app.ImprovedBiasLabel.Text = 'After Correction: --';
            app.ImprovedSDLabel.Text   = '';
            app.CorrQualityLabel.Text  = '';

            slopeDev = abs(app.Stats.slope - 1.0);
            if slopeDev > 0.5
                app.SlopeWarningLabel.Text = sprintf('⚠ Slope = %.2f: range-dependent disagreement detected; proportional models may be informative.', app.Stats.slope);
                app.SlopeWarningLabel.FontColor = [0.75 0.45 0.0];
                app.SlopeWarningLabel.Visible = 'on';
            elseif slopeDev > 0.25
                app.SlopeWarningLabel.Text = sprintf('⚠ Slope = %.2f: possible proportional disagreement.', app.Stats.slope);
                app.SlopeWarningLabel.FontColor = [0.75 0.45 0.0];
                app.SlopeWarningLabel.Visible = 'on';
            else
                app.SlopeWarningLabel.Visible = 'off';
            end
            app.CorrectionModel        = struct();  
            app.SmallNWarningLabel.Visible = 'off';
            app.CorrQualityLabel.Text  = '';
            
            app.ShowCorrectedSwitch.Value = false;
            app.ShowCorrectedSwitch.Enable = 'off';
            app.ShowOriginalCDISwitch.Value = true;
            app.KeepCorrectionSwitch.Value = false;
            app.KeepCorrectionSwitch.Enable = 'off';
            
            app.CorrViewSwitch.Value = false;
            app.CorrViewSwitch.Enable = 'off';
            app.BAViewSwitch.Value = false;
            app.BAViewSwitch.Enable = 'off';
            
            app.ComparePlotsButton.Enable = 'off'; 
            
            if ~isempty(app.CorrectedTrendLine) && isvalid(app.CorrectedTrendLine)
                delete(app.CorrectedTrendLine);
                app.CorrectedTrendLine = [];
            end
            if ~isempty(app.KeptCorrectionLines)
                for ki = 1:numel(app.KeptCorrectionLines)
                    if isvalid(app.KeptCorrectionLines{ki})
                        delete(app.KeptCorrectionLines{ki});
                    end
                end
                app.KeptCorrectionLines = {};
            end
            app.KeptCorrectionData = {};

            cla(app.TimeAxes);
            hold(app.TimeAxes, 'on');
            
            raw_abl_t = raw_abl_t_pre;
            raw_abl_y = raw_abl_y_pre;
            valid_abl = valid_abl_pre;

            full_cdi_t = full_cdi_t_pre;
            full_cdi_y = full_cdi_y_pre;
            valid_cdi  = valid_cdi_pre;

            if hasOverlap
                pad = minutes(max(30, minutes(overlap_end - overlap_start) * 0.02));
                xlim_start = overlap_start - pad;
                xlim_end   = overlap_end   + pad;
            else
                pad = minutes(30);
                xlim_start = min(abl_start, cdi_start);
                xlim_end   = max(abl_end,   cdi_end);
            end

            abl_in_window  = valid_abl & raw_abl_t >= cdi_start & raw_abl_t <= cdi_end;
            abl_out_window = valid_abl & ~abl_in_window;
            
            cdi_in_window = valid_cdi & full_cdi_t >= overlap_start - pad & ...
                            full_cdi_t <= overlap_end + pad;

            p2 = plot(app.TimeAxes, full_cdi_t(cdi_in_window), full_cdi_y(cdi_in_window), '-', ...
                'DisplayName', 'CDI Continuous', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
            app.OriginalCDILine = p2;

            if any(abl_in_window)
                p1 = plot(app.TimeAxes, raw_abl_t(abl_in_window), raw_abl_y(abl_in_window), 'o', ...
                    'DisplayName', 'ABL (ref draws)', 'LineWidth', 1.5, 'MarkerSize', 7, ...
                    'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerFaceColor', [0 0.4470 0.7410]);
                p1.DataTipTemplate.DataTipRows(1).Label = 'Time';
                p1.DataTipTemplate.DataTipRows(2).Label = 'ABL';
            end

            if any(abl_out_window)
                p0 = plot(app.TimeAxes, raw_abl_t(abl_out_window), raw_abl_y(abl_out_window), 'o', ...
                    'DisplayName', 'ABL (no CDI coverage)', 'LineWidth', 1, 'MarkerSize', 6, ...
                    'MarkerEdgeColor', [0.6 0.6 0.6], 'MarkerFaceColor', 'none', ...
                    'LineStyle', 'none');
                p0.DataTipTemplate.DataTipRows(1).Label = 'Time';
                p0.DataTipTemplate.DataTipRows(2).Label = 'ABL (unpaired)';
            end

            p2.DataTipTemplate.DataTipRows(1).Label = 'Time';
            p2.DataTipTemplate.DataTipRows(2).Label = 'CDI';

            if hasOverlap
                xlim(app.TimeAxes, [xlim_start, xlim_end]);
            end

            if app.FitWindowCheckBox.Value && ~isempty(strtrim(app.FitWindowStartEdit.Value))
                try
                    winS = datetime(app.FitWindowStartEdit.Value, 'InputFormat', 'dd.MM.yyyy HH:mm');
                    winE = datetime(app.FitWindowEndEdit.Value,   'InputFormat', 'dd.MM.yyyy HH:mm');
                    cdi_win_vals = full_cdi_y(valid_cdi & full_cdi_t >= winS & full_cdi_t <= winE);
                    abl_win_vals = raw_abl_y(valid_abl & raw_abl_t >= winS & raw_abl_t <= winE);
                    all_win = [cdi_win_vals; abl_win_vals];
                    all_win = all_win(~isnan(all_win));
                    if numel(all_win) >= 2
                        y_lo = min(all_win);
                        y_hi = max(all_win);
                        y_span = max(y_hi - y_lo, 1);
                        ylim(app.TimeAxes, [y_lo - y_span*0.15, y_hi + y_span*0.15]);
                    end
                catch
                end
            end

            patLbl = app.PatientIDDropDown.Value;
            if ~hasOverlap
                title(app.TimeAxes, [param ' — ' patLbl ' ⚠ No ABL-CDI overlap detected']);
            else
                title(app.TimeAxes, [param ' — Patient: ' patLbl]);
            end

            legend(app.TimeAxes, 'Location', 'best');
            xlabel(app.TimeAxes, 'Time');
            ylabel(app.TimeAxes, param);
            grid(app.TimeAxes, 'on');
            hold(app.TimeAxes, 'off');

            if hasOverlap && (~app.FitWindowCheckBox.Value || ...
                    isempty(strtrim(app.FitWindowStartEdit.Value)))
                app.FitWindowStartEdit.Value = datestr(overlap_start, 'dd.mm.yyyy HH:MM');
                app.FitWindowEndEdit.Value   = datestr(overlap_end,   'dd.mm.yyyy HH:MM');
            end

            drawFitWindowShade(app);

            CorrViewSwitchChanged(app, []);
            BAViewSwitchChanged(app, []);

            drawnow;  

            app.ExportButton.Enable = 'on';
            app.ExportFigureButton.Enable = 'on';
            app.ApplyCorrectionButton.Enable = 'on';
            app.StatusLabel.Text = sprintf('Analysis complete. %d pairs matched.', app.Stats.N);
        end

        function ExportButtonPushed(app, ~)
            if isempty(app.Stats) || ~isfield(app.Stats, 'N')
                uialert(app.UIFigure, 'No analysis to export. Please run analysis first.', 'No Data');
                return;
            end
            [file, path] = uiputfile({'*.xlsx';'*.csv'}, 'Save Results As');
            if isequal(file,0), return; end
            fullpath = fullfile(path, file);
            try
                statsTable = struct2table(app.Stats);
                writetable(statsTable, fullpath, 'Sheet', 'Statistics');
                if ~isempty(app.Aligned_Data)
                    dataTable = timetable2table(app.Aligned_Data);
                    writetable(dataTable, fullpath, 'Sheet', 'PairedData');
                end
                app.StatusLabel.Text = ['Exported to: ' file];
            catch ME
                uialert(app.UIFigure, ['Export failed: ' ME.message], 'Export Error');
            end
        end
        
        function ExportFigureButtonPushed(app, ~)
            if isempty(app.Stats) || ~isfield(app.Stats, 'N')
                uialert(app.UIFigure, 'No analysis to export. Please run analysis first.', 'No Data');
                return;
            end

            param   = app.CurrentParam;
            patID   = app.PatientIDDropDown.Value;
            safePat = regexprep(string(patID), '[^a-zA-Z0-9_]', '_');
            safeParam = regexprep(param, '[^a-zA-Z0-9_]', '_');
            defaultName = sprintf('%s_%s_figures', safePat, safeParam);

            [file, path] = uiputfile({'*.svg', 'SVG vector (*.svg)'}, ...
                'Save Figure As', defaultName);
            if isequal(file, 0), return; end

            [~, baseName, ~] = fileparts(file);

            try
                param   = app.CurrentParam;
                xABL    = app.Stats.xABL;
                yOrig   = app.Stats.yCDI;
                biasO   = app.Stats.bias;
                sdO     = app.Stats.sd;
                hasCorrected = isfield(app.CorrectionModel,'yCorrected') && ~isempty(app.CorrectionModel.yCorrected);
                if hasCorrected
                    yCorr = app.CorrectionModel.yCorrected;
                    if numel(yCorr)~=numel(xABL)
                        nm=min(numel(yCorr),numel(xABL));
                        yCorr=yCorr(1:nm); xABL=xABL(1:nm); yOrig=yOrig(1:nm);
                    end
                    
                    % Filter baseline so Before/After compare the exact same pairs
                    validPairs = ~isnan(yCorr) & ~isnan(xABL);
                    xABL = xABL(validPairs);
                    yOrig = yOrig(validPairs);
                    yCorr = yCorr(validPairs);
                    
                    biasO = mean(yOrig - xABL, 'omitnan');
                    sdO = std(yOrig - xABL, 'omitnan');
                    
                    diffCorr = yCorr-xABL;
                    biasC = mean(diffCorr,'omitnan'); sdC = std(diffCorr,'omitnan');
                    loaUpC = biasC+1.96*sdC; loaLoC = biasC-1.96*sdC;
                    methodLbl = app.CorrectionMethodDropDown.Value;
                    if isfield(app.CorrectionModel,'autoSelected') && app.CorrectionModel.autoSelected
                        methodLbl = app.CorrectionModel.autoWinner;
                    end
                end

                allVals = [xABL; yOrig];
                if hasCorrected, allVals=[allVals; yCorr]; end
                vv=allVals(~isnan(allVals)&~isinf(allVals));
                if isempty(vv), sMin=0; sMax=1;
                else
                    sp=max(vv)-min(vv); if sp==0||isnan(sp),sp=0.1; end
                    sMin=min(vv)-sp*0.1; sMax=max(vv)+sp*0.1;
                end

                diffO = yOrig-xABL;
                avgO  = (xABL+yOrig)/2;
                allAvg=avgO; allDiff=diffO;
                if hasCorrected
                    allAvg=[allAvg;(xABL+yCorr)/2]; allDiff=[allDiff;diffCorr];
                end
                va=allAvg(~isnan(allAvg)&~isinf(allAvg));
                if isempty(va), bxl=[0 1];
                else, sp=max(va)-min(va); if sp==0,sp=0.1; end
                    bxl=[min(va)-sp*0.15, max(va)+sp*0.15];
                end
                vd=allDiff(~isnan(allDiff)&~isinf(allDiff));
                if isempty(vd), byl=[-1 1];
                else, sp=max(vd)-min(vd); if sp==0,sp=0.1; end
                    byl=[min(vd)-sp*0.25, max(vd)+sp*0.25];
                end
                loaUpO=biasO+1.96*sdO; loaLoO=biasO-1.96*sdO;

                hFig = figure('Visible','off','Units','centimeters','Position',[0 0 44 26]);
                set(hFig,'Color','white');

                axT = copyobj(app.TimeAxes, hFig);
                axT.Units='normalized';
                axT.Position=[0.05 0.55 0.92 0.40];
                axT.FontSize=11;
                axT.XLabel.FontSize=11;
                axT.YLabel.FontSize=11;
                axT.Title.FontSize=12;
                legend(axT,'Location','northeast','FontSize',9);

                axCB = axes(hFig,'Position',[0.05 0.08 0.18 0.36]);
                hold(axCB,'on');
                scatter(axCB,xABL,yOrig,30,'filled','MarkerFaceColor',[0 0.4470 0.7410],...
                    'DisplayName',sprintf('r=%.3f',app.Stats.r));
                plot(axCB,[sMin sMax],[sMin sMax],'k--','LineWidth',0.8,'DisplayName','Identity');
                if ~any(isnan(app.Stats.p))
                    xf=linspace(sMin,sMax,100);
                    plot(axCB,xf,polyval(app.Stats.p,xf),'-','Color',[0 0.4470 0.7410],'LineWidth',1.2,'DisplayName','Regression');
                end
                axis(axCB,'square'); xlim(axCB,[sMin sMax]); ylim(axCB,[sMin sMax]);
                xlabel(axCB,['ABL ' param],'FontSize',8); ylabel(axCB,['CDI ' param],'FontSize',8);
                title(axCB,sprintf('Correlation Before\nr=%.3f  N=%d',app.Stats.r,numel(xABL)),'FontSize',8);
                legend(axCB,'Location','southeast','FontSize',6); grid(axCB,'on'); hold(axCB,'off');

                axCA = axes(hFig,'Position',[0.28 0.08 0.18 0.36]);
                hold(axCA,'on');
                if hasCorrected
                    scatter(axCA,xABL,yCorr,30,'filled','MarkerFaceColor',[0.2 0.7 0.2],...
                        'DisplayName',sprintf('r=%.3f',app.CorrectionModel.r_new));
                    plot(axCA,[sMin sMax],[sMin sMax],'k--','LineWidth',0.8,'DisplayName','Identity');
                    if ~any(isnan(app.CorrectionModel.pCorr))
                        xf=linspace(sMin,sMax,100);
                        plot(axCA,xf,polyval(app.CorrectionModel.pCorr,xf),'-','Color',[0.1 0.5 0.1],'LineWidth',1.2,'DisplayName','Regression');
                    end
                    shortLbl = methodLbl;
                    if numel(shortLbl) > 18, shortLbl = [shortLbl(1:16) '...']; end
                    title(axCA,sprintf('Correlation After\n(%s)  r=%.3f',shortLbl,app.CorrectionModel.r_new),'FontSize',8);
                    legend(axCA,'Location','southeast','FontSize',6);
                else
                    text(axCA,0.5,0.5,'No correction applied','Units','normalized',...
                        'HorizontalAlignment','center','FontSize',8,'Color',[0.5 0.5 0.5]);
                    title(axCA,'Correlation After','FontSize',8);
                end
                axis(axCA,'square'); xlim(axCA,[sMin sMax]); ylim(axCA,[sMin sMax]);
                xlabel(axCA,['ABL ' param],'FontSize',8); ylabel(axCA,['CDI ' param],'FontSize',8);
                grid(axCA,'on'); hold(axCA,'off');

                axBB = axes(hFig,'Position',[0.53 0.08 0.20 0.36]);
                hold(axBB,'on');
                scatter(axBB,avgO,diffO,30,'filled','MarkerFaceColor',[0 0.4470 0.7410],...
                    'DisplayName',sprintf('N=%d',numel(diffO)));
                plot(axBB,bxl,[0 0],'k:','LineWidth',0.8,'HandleVisibility','off');
                plot(axBB,bxl,[biasO biasO],'-','Color',[0 0.4470 0.7410],'LineWidth',1.5,...
                    'DisplayName',sprintf('Bias=%.4f',biasO));
                plot(axBB,bxl,[loaUpO loaUpO],'--','Color',[0 0.4470 0.7410],'LineWidth',1,...
                    'DisplayName',sprintf('LoA=[%.3f,%.3f]',loaLoO,loaUpO));
                plot(axBB,bxl,[loaLoO loaLoO],'--','Color',[0 0.4470 0.7410],'LineWidth',1,'HandleVisibility','off');
                xlim(axBB,bxl); ylim(axBB,byl);
                xlabel(axBB,['Mean  ' param],'FontSize',8);
                ylabel(axBB,'CDI - ABL','FontSize',8);
                title(axBB,sprintf('Bland-Altman Before\nBias=%.4f  SD=%.4f',biasO,sdO),'FontSize',8);
                legend(axBB,'Location','northeast','FontSize',6); grid(axBB,'on'); hold(axBB,'off');

                axBA = axes(hFig,'Position',[0.77 0.08 0.20 0.36]);
                hold(axBA,'on');
                if hasCorrected
                    avgC=(xABL+yCorr)/2;
                    scatter(axBA,avgC,diffCorr,30,'filled','MarkerFaceColor',[0.2 0.7 0.2],...
                        'DisplayName',sprintf('N=%d',numel(diffCorr)));
                    plot(axBA,bxl,[0 0],'k:','LineWidth',0.8,'HandleVisibility','off');
                    plot(axBA,bxl,[biasC biasC],'-','Color',[0.1 0.5 0.1],'LineWidth',1.5,...
                        'DisplayName',sprintf('Bias=%.4f',biasC));
                    plot(axBA,bxl,[loaUpC loaUpC],'--','Color',[0.1 0.5 0.1],'LineWidth',1,...
                        'DisplayName',sprintf('LoA=[%.3f,%.3f]',loaLoC,loaUpC));
                    plot(axBA,bxl,[loaLoC loaLoC],'--','Color',[0.1 0.5 0.1],'LineWidth',1,'HandleVisibility','off');
                    shortLbl2 = methodLbl;
                    if numel(shortLbl2) > 18, shortLbl2 = [shortLbl2(1:16) '...']; end
                    title(axBA,sprintf('Bland-Altman After\n(%s)',shortLbl2),'FontSize',8);
                    legend(axBA,'Location','northeast','FontSize',6);
                else
                    text(axBA,0.5,0.5,'No correction applied','Units','normalized',...
                        'HorizontalAlignment','center','FontSize',8,'Color',[0.5 0.5 0.5]);
                    title(axBA,'Bland-Altman After','FontSize',8);
                end
                xlim(axBA,bxl); ylim(axBA,byl);
                xlabel(axBA,['Mean  ' param],'FontSize',8);
                ylabel(axBA,'CDI - ABL','FontSize',8);
                grid(axBA,'on'); hold(axBA,'off');

                outPathAll = fullfile(path,[baseName '_report.svg']);
                print(hFig, outPathAll, '-dsvg');
                close(hFig);

                hF2 = figure('Visible','off','Units','centimeters','Position',[0 0 30 9]);
                set(hF2,'Color','white');
                axT2 = copyobj(app.TimeAxes, hF2);
                axT2.Units='normalized'; axT2.Position=[0.07 0.15 0.90 0.76];
                axT2.FontSize=12;
                axT2.XLabel.FontSize=12;
                axT2.YLabel.FontSize=12;
                axT2.Title.FontSize=13;
                axT2.Title.String = sprintf('%s  -  Patient: %s', param, app.PatientIDDropDown.Value);
                legend(axT2,'Location','northeast','FontSize',10);
                print(hF2, fullfile(path,[baseName '_timeseries.svg']), '-dsvg');
                close(hF2);

                hF3=figure('Visible','off','Units','centimeters','Position',[0 0 12 12],'Color','white');
                ax3=copyobj(app.CorrelationAxes,hF3); ax3.Units='normalized'; ax3.Position=[0.14 0.14 0.78 0.78];
                print(hF3, fullfile(path,[baseName '_correlation.svg']), '-dsvg');
                close(hF3);

                hF4=figure('Visible','off','Units','centimeters','Position',[0 0 12 10],'Color','white');
                ax4=copyobj(app.BlandAltmanAxes,hF4); ax4.Units='normalized'; ax4.Position=[0.14 0.14 0.78 0.78];
                print(hF4, fullfile(path,[baseName '_blandaltman.svg']), '-dsvg');
                close(hF4);

                hF5 = figure('Visible','off','Units','centimeters','Position',[0 0 8 7]);
                set(hF5,'Color','white');
                ax5 = axes(hF5,'Visible','off');
                ax5.Units='normalized'; ax5.Position=[0 0 1 1];

                statsLines = {};
                statsLines{end+1} = sprintf('N Pairs: %d', app.Stats.N);
                statsLines{end+1} = sprintf('Bias:    %.4f', app.Stats.bias);
                statsLines{end+1} = sprintf('SD:      %.4f', app.Stats.sd);
                statsLines{end+1} = sprintf('95%% LoA: [%.4f, %.4f]', app.Stats.loa_lo, app.Stats.loa_up);
                statsLines{end+1} = sprintf('r =      %.4f', app.Stats.r);
                statsLines{end+1} = sprintf('y = %.3fx %+.3f', app.Stats.slope, app.Stats.intercept);

                if isfield(app.CorrectionModel,'yCorrected') && ~isempty(app.CorrectionModel.yCorrected)
                    xABLs = app.Stats.xABL;
                    yCorrs = app.CorrectionModel.yCorrected;
                    if numel(yCorrs) ~= numel(xABLs)
                        nm = min(numel(yCorrs),numel(xABLs));
                        yCorrs = yCorrs(1:nm); xABLs = xABLs(1:nm);
                    end
                    dC = yCorrs - xABLs;
                    biasC2 = mean(dC,'omitnan'); sdC2 = std(dC,'omitnan');
                    loaUpC2 = biasC2+1.96*sdC2; loaLoC2 = biasC2-1.96*sdC2;
                    mLbl = app.CorrectionMethodDropDown.Value;
                    if isfield(app.CorrectionModel,'autoSelected') && app.CorrectionModel.autoSelected
                        mLbl = app.CorrectionModel.autoWinner;
                    end
                    statsLines{end+1} = '';
                    statsLines{end+1} = sprintf('After %s:', mLbl);
                    statsLines{end+1} = sprintf('Bias:    %.4f', biasC2);
                    statsLines{end+1} = sprintf('SD:      %.4f', sdC2);
                    statsLines{end+1} = sprintf('95%% LoA: [%.4f, %.4f]', loaLoC2, loaUpC2);
                    if isfield(app.CorrectionModel,'r_new')
                        statsLines{end+1} = sprintf('r =      %.4f', app.CorrectionModel.r_new);
                    end
                end

                yStep = 0.80 / (numel(statsLines) + 1);
                for si = 1:numel(statsLines)
                    yPos = 0.84 - (si-1)*yStep;
                    isHeader = ~isempty(statsLines{si}) && statsLines{si}(end)==':';
                    fw = 'normal'; if isHeader, fw='bold'; end
                    text(ax5, 0.05, yPos, statsLines{si}, ...
                        'Units','normalized','FontSize',9,'FontName','Consolas',...
                        'FontWeight',fw,'Color',[0.1 0.1 0.1],'Interpreter','none');
                end
                rectangle('Parent',ax5,'Position',[0.01 0.01 0.98 0.98],...
                    'EdgeColor',[0.4 0.4 0.4],'LineWidth',1);
                text(ax5, 0.05, 0.96, sprintf('Statistics - %s', app.CurrentParam), ...
                    'Units','normalized','FontSize',9,'FontWeight','bold',...
                    'Color',[0.1 0.1 0.5],'Interpreter','none','FontName','Consolas');
                text(ax5, 0.05, 0.90, sprintf('Patient: %s', app.PatientIDDropDown.Value), ...
                    'Units','normalized','FontSize',8,'FontWeight','normal',...
                    'Color',[0.2 0.2 0.6],'Interpreter','none','FontName','Consolas');

                outPathStats = fullfile(path,[baseName '_statistics.svg']);
                print(hF5, outPathStats, '-dsvg');
                close(hF5);

                app.StatusLabel.Text = sprintf('Figures saved: %s', baseName);
                uialert(app.UIFigure, ...
                    sprintf('Saved to:\n%s\n\n  %s_report.svg  (combined)\n  %s_timeseries.svg\n  %s_correlation.svg\n  %s_blandaltman.svg\n  %s_statistics.svg', ...
                    path, baseName, baseName, baseName, baseName, baseName), ...
                    'Export Complete', 'Icon', 'success');

            catch ME
                try; close(hFig);  catch; end
                try; close(hF2);   catch; end
                try; close(hF3);   catch; end
                try; close(hF4);   catch; end
                try; close(hF5);   catch; end
                uialert(app.UIFigure, ['Figure export failed: ' ME.message], 'Export Error');
            end
        end

        % --- MODEL APPLICATION / 7-MODEL LOO-CV SELECTION ---
        function ApplyCorrectionButtonPushed(app, ~)
            if isempty(app.Aligned_Data) || ~isfield(app.Stats, 'TimeShift')
                uialert(app.UIFigure, 'Please run analysis first.', 'No Data');
                return;
            end
            
            method = app.CorrectionMethodDropDown.Value;
            param  = app.CurrentParam;
            
            if ~isempty(app.CorrectedTrendLine) && isvalid(app.CorrectedTrendLine)
                delete(app.CorrectedTrendLine);
                app.CorrectedTrendLine = [];
            end
            app.ShowCorrectedSwitch.Value = false;
            
            varNames = app.Aligned_Data.Properties.VariableNames;
            ablCol = ''; cdiCol = '';
            for i = 1:length(varNames)
                v = varNames{i};
                if contains(v, param) && contains(v, 'abl')
                    ablCol = v;
                end
                if contains(v, param) && contains(v, 'cdi')
                    cdiCol = v;
                end
            end
            
            if isempty(ablCol) || isempty(cdiCol)
                uialert(app.UIFigure, 'Could not find data columns.', 'Error');
                return;
            end
            
            xABL = app.Aligned_Data.(ablCol);
            yCDI = app.Aligned_Data.(cdiCol);
            validIdx = ~isnan(xABL) & ~isnan(yCDI);
            xABL = xABL(validIdx);
            yCDI = yCDI(validIdx);
            pairedTimes = app.Aligned_Data.Time(validIdx);

            fitWinMask = getFitWindowMask(app, pairedTimes);
            xABL_fit = xABL(fitWinMask);
            yCDI_fit = yCDI(fitWinMask);

            app.CorrectionModel.autoSelected = false;
            app.CorrectionModel.fitWindowUsed = app.FitWindowCheckBox.Value;
            app.CorrectionModel.nFitPairs = numel(xABL_fit);
            yCorrected = []; 

            shiftMins   = app.Stats.TimeShift;
            fullCDITime = app.CDI_Table.Time + minutes(shiftMins);
            fullCDIVals = app.CDI_Table.(param);
            
            switch method
                case 'Hybrid (Time-Series + Deming)'
                    tau_r = app.TauSpinner.Value;
                    tau_f = app.TauFallSpinner.Value;
                    lam   = app.DemingLambdaEditField.Value;
                    w1    = app.SmoothW1Spinner.Value;
                    
                    cdi_fast_full = computeAsymmetricFastCDI(app, fullCDIVals, fullCDITime, w1, tau_r, tau_f);
                    validCDIMask  = ~isnan(fullCDIVals);
                    searchTimes   = fullCDITime(validCDIMask);
                    searchVals    = cdi_fast_full(validCDIMask);
                    cdi_fast_paired = zeros(size(xABL));
                    for k = 1:numel(xABL)
                        [~, bestIdx] = min(abs(searchTimes - pairedTimes(k)));
                        cdi_fast_paired(k) = searchVals(bestIdx);
                    end
                    cdi_fast_paired_fit = cdi_fast_paired(fitWinMask);
                    cleanMask = robustCleanMask(app, xABL_fit, cdi_fast_paired_fit);
                    xClean = xABL_fit(cleanMask);
                    yClean = cdi_fast_paired_fit(cleanMask);
                    
                    [slope, intercept] = fitWeightedDeming(app, xClean, yClean, lam);
                    yCorrected = (cdi_fast_paired(fitWinMask) - intercept) / slope;
                    yCorrected(~cleanMask) = NaN;
                    
                    app.CorrectionModel.type = 'hybrid';
                    app.CorrectionModel.tau_rise = tau_r;
                    app.CorrectionModel.tau_fall = tau_f;
                    app.CorrectionModel.lam = lam;
                    app.CorrectionModel.slope = slope;
                    app.CorrectionModel.intercept = intercept;
                    app.CorrectionModel.formula = sprintf('Step 1: CDI_s = movmean(CDI_raw, [%d 0])\nStep 2: dCDI/dt = movmean(diff(CDI_s)/dt, [%d 0])\nStep 3: CDI_fast = CDI_s + [τ_r=%.1f / τ_f=%.1f] * (dCDI/dt)\nStep 4: CDI_fast = movmean(CDI_fast, [2 0])\nStep 5: CDI_corrected = (CDI_fast - %.4f) / %.4f  [Linnet Weighted Deming λ=%.2f]\nHybrid filtering is causal: current output uses current and prior CDI samples only.', max(0, w1-1), max(0, w1-1), tau_r, tau_f, intercept, slope, lam);

                case 'Weighted Deming (Linnet)'
                    cleanMask = robustCleanMask(app, xABL_fit, yCDI_fit);
                    lam = app.DemingLambdaEditField.Value;
                    xC = xABL_fit(cleanMask); yC = yCDI_fit(cleanMask);
                    [slope, intercept] = fitWeightedDeming(app, xC, yC, lam);
                    yCorrected = (yCDI_fit - intercept) / slope;
                    yCorrected(~cleanMask) = NaN;
                    app.CorrectionModel.type = 'weighted_deming';
                    app.CorrectionModel.slope = slope;
                    app.CorrectionModel.intercept = intercept;
                    app.CorrectionModel.lam = lam;
                    app.DemingLambdaEditField.Value = lam;
                    app.CorrectionModel.formula = sprintf('Model: Raw_CDI = %.4f*ABL %+.4f  [Linnet Weighted Deming λ=%.2f]\nCorrected = (Raw_CDI %+.4f) / %.4f\n(Fitted on %d/%d window pairs, %d robust [|diff - med| <= 4.5*MAD])', slope, intercept, lam, -intercept, slope, numel(xABL_fit), numel(xABL), sum(cleanMask));

                case 'Bias Correction'
                    cleanMask = robustCleanMask(app, xABL_fit, yCDI_fit);
                    bias = mean(yCDI_fit(cleanMask) - xABL_fit(cleanMask), 'omitnan');
                    yCorrected = yCDI_fit - bias;
                    yCorrected(~cleanMask) = NaN;
                    app.CorrectionModel.type = 'bias';
                    app.CorrectionModel.bias = bias;
                    app.CorrectionModel.formula = sprintf('Model: Raw_CDI = ABL %+.4f\nCorrected = Raw_CDI %+.4f\n(Fitted on %d/%d window pairs, %d robust [|diff - med| <= 4.5*MAD])', bias, -bias, numel(xABL_fit), numel(xABL), sum(cleanMask));
                    
                case 'OLS (ABL is X)'
                    cleanMask = robustCleanMask(app, xABL_fit, yCDI_fit);
                    xC = xABL_fit(cleanMask); yC = yCDI_fit(cleanMask);
                    p_fit = polyfit(xC, yC, 1);
                    slope = p_fit(1); intercept = p_fit(2);
                    if abs(slope) < 1e-4
                        slope = 1; intercept = mean(yC) - mean(xC);
                    end
                    yCorrected = (yCDI_fit - intercept) / slope;
                    yCorrected(~cleanMask) = NaN;
                    app.CorrectionModel.type = 'ols_abl_x';
                    app.CorrectionModel.slope = slope;
                    app.CorrectionModel.intercept = intercept;
                    app.CorrectionModel.formula = sprintf('Model: Raw_CDI = %.4f*ABL %+.4f\nCorrected = (Raw_CDI %+.4f) / %.4f\n(Fitted on %d/%d window pairs, %d robust [|diff - med| <= 4.5*MAD])', slope, intercept, -intercept, slope, numel(xABL_fit), numel(xABL), sum(cleanMask));

                case 'Proportional Correction'
                    cleanMask = robustCleanMask(app, xABL_fit, yCDI_fit);
                    xC = xABL_fit(cleanMask); yC = yCDI_fit(cleanMask);
                    ratio = mean(xC ./ yC, 'omitnan');
                    yCorrected = yCDI_fit * ratio;
                    yCorrected(~cleanMask) = NaN;
                    app.CorrectionModel.type = 'proportional';
                    app.CorrectionModel.ratio = ratio;
                    app.CorrectionModel.formula = sprintf('Ratio (ABL/CDI) = %.4f\nCorrected = Raw_CDI * %.4f\n(Fitted on %d/%d window pairs, %d robust [|diff - med| <= 4.5*MAD])', ratio, ratio, numel(xABL_fit), numel(xABL), sum(cleanMask));
                    
                case 'Deming Regression'
                    cleanMask = robustCleanMask(app, xABL_fit, yCDI_fit);
                    xC = xABL_fit(cleanMask); yC = yCDI_fit(cleanMask);
                    lam = app.DemingLambdaEditField.Value;
                    n2 = numel(xC);
                    xm = mean(xC,'omitnan'); ym = mean(yC,'omitnan');
                    sxx = sum((xC-xm).^2,'omitnan')/(n2-1);
                    syy = sum((yC-ym).^2,'omitnan')/(n2-1);
                    sxy = sum((xC-xm).*(yC-ym),'omitnan')/(n2-1);
                    denom = 2*sxy;
                    if abs(denom)<1e-10, slope=1;
                    else, slope=(syy-lam*sxx+sqrt((syy-lam*sxx)^2+4*lam*sxy^2))/denom;
                    end
                    if isnan(slope)||isinf(slope)||abs(slope)<1e-4
                        slope=1; intercept=ym-xm;
                    else, intercept=ym-slope*xm;
                    end
                    yCorrected = (yCDI_fit - intercept) / slope;
                    yCorrected(~cleanMask) = NaN;
                    app.CorrectionModel.type = 'deming';
                    app.CorrectionModel.slope = slope;
                    app.CorrectionModel.intercept = intercept;
                    app.CorrectionModel.lam = lam;
                    app.CorrectionModel.formula = sprintf('Model: Raw_CDI = %.4f*ABL %+.4f  [Standard Deming λ=%.2f]\nCorrected = (Raw_CDI %+.4f) / %.4f\n(Fitted on %d/%d window pairs, %d robust [|diff - med| <= 4.5*MAD])', slope, intercept, lam, -intercept, slope, numel(xABL_fit), numel(xABL), sum(cleanMask));

                case 'Passing-Bablok'
                    cleanMask = robustCleanMask(app, xABL_fit, yCDI_fit);
                    xC = xABL_fit(cleanMask); yC = yCDI_fit(cleanMask);
                    n2 = numel(xC);
                    slopes_pb = zeros(n2*(n2-1)/2, 1);
                    idx_pb = 1;
                    for i2 = 1:n2-1
                        for j2 = i2+1:n2
                            if xC(j2) ~= xC(i2)
                                slopes_pb(idx_pb) = (yC(j2)-yC(i2))/(xC(j2)-xC(i2));
                            else, slopes_pb(idx_pb) = NaN;
                            end
                            idx_pb = idx_pb + 1;
                        end
                    end
                    slopes_pb = slopes_pb(~isnan(slopes_pb));
                    slope = median(slopes_pb,'omitnan');
                    if isnan(slope)||isinf(slope)||abs(slope)<1e-4, slope=1; end
                    intercept = median(yC - slope.*xC,'omitnan');
                    if isnan(intercept), intercept=mean(yC)-mean(xC); end
                    yCorrected = (yCDI_fit - intercept) / slope;
                    yCorrected(~cleanMask) = NaN;
                    app.CorrectionModel.type = 'passing-bablok';
                    app.CorrectionModel.slope = slope;
                    app.CorrectionModel.intercept = intercept;
                    app.CorrectionModel.formula = sprintf('Model: Raw_CDI = %.4f*ABL %+.4f  [Simplified Passing-Bablok]\nCorrected = (Raw_CDI %+.4f) / %.4f\n(Fitted on %d/%d window pairs, %d robust [|diff - med| <= 4.5*MAD])', slope, intercept, -intercept, slope, numel(xABL_fit), numel(xABL), sum(cleanMask));

                case 'Auto (Best Model)'
                    app.StatusLabel.Text = 'Auto: running strict out-of-sample LOO-CV (7 candidate models)...'; drawnow;

                    % Corrected literal LaTeX slashes to Unicode so UI doesn't break rendering
                    candidateNames   = {'Bias Correction', ...
                                        'OLS (ABL is X)', ...
                                        'Proportional Correction', ...
                                        'Deming Regression (fixed λ=1)', ...
                                        'Weighted Deming (Linnet, tuned λ)', ...
                                        'Passing-Bablok (simplified)', ...
                                        'Hybrid (Time-Series + Deming)'};
                    candidateRMSE    = nan(1, 7);
                    candidateLoASpan = nan(1, 7);
                    n = numel(xABL_fit);

                    lam_grid  = [0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0];
                    tau_grid  = 0:1.0:8;
                    w1_grid   = [4, 8, 12, 16, 24, 32];

                    validCDIMask = ~isnan(fullCDIVals);
                    searchTimes  = fullCDITime(validCDIMask);
                    fitTimes     = pairedTimes(fitWinMask);

                    if numel(fitTimes) >= 2
                        winStartH = min(fitTimes);
                        winEndH   = max(fitTimes);
                        winMaskCDIH = searchTimes >= winStartH & searchTimes <= winEndH;
                    else
                        winMaskCDIH = true(size(searchTimes));
                    end

                    % Pre-compute fast dynamic CDI traces across candidate grids to accelerate LOO folds
                    nTau = numel(tau_grid); 
                    nW1  = numel(w1_grid);
                    fastPairedFit = zeros(n, nTau, nW1);
                    fastRoughnessWin = zeros(nTau, nW1);

                    for wi = 1:nW1
                        w1v = w1_grid(wi);
                        for ti = 1:nTau
                            tv = tau_grid(ti);
                            cdi_fast_full = computeAsymmetricFastCDI(app, fullCDIVals, fullCDITime, w1v, tv, tv);
                            sv = cdi_fast_full(validCDIMask);
                            
                            for k = 1:n
                                [~, bi] = min(abs(searchTimes - fitTimes(k)));
                                fastPairedFit(k, ti, wi) = sv(bi);
                            end
                            
                            sv_win = sv(winMaskCDIH);
                            rv = sv_win(~isnan(sv_win));
                            if numel(rv) > 2
                                fastRoughnessWin(ti, wi) = std(diff(rv), 'omitnan');
                            else
                                fastRoughnessWin(ti, wi) = 1.0;
                            end
                        end
                    end

                    raw_valid_h = fullCDIVals(validCDIMask);
                    raw_valid_h_win = raw_valid_h(winMaskCDIH);
                    if numel(raw_valid_h_win) > 2
                        raw_rough = std(diff(raw_valid_h_win), 'omitnan');
                        if raw_rough == 0 || isnan(raw_rough), raw_rough = 1.0; end
                    else
                        raw_rough = 1.0;
                    end

                    % =========================================================================
                    % STRICT LEAVE-ONE-OUT CROSS-VALIDATION (Hyperparameter tuning strictly in-fold)
                    % =========================================================================
                    for ci = 1:7
                        looErrors = nan(n, 1);
                        for i = 1:n
                            leaveIdx = [1:i-1, i+1:n];
                            xTr = xABL_fit(leaveIdx);
                            yTr = yCDI_fit(leaveIdx);
                            xTe = xABL_fit(i);
                            yTe = yCDI_fit(i);

                            try
                                cMask = robustCleanMask(app, xTr, yTr);
                                xTrC  = xTr(cMask); 
                                yTrC  = yTr(cMask);

                                switch ci
                                    case 1 % Bias Correction
                                        b = mean(yTrC - xTrC, 'omitnan');
                                        pred = yTe - b;

                                    case 2 % OLS Regression
                                        if numel(xTrC) >= 2 && std(xTrC, 'omitnan') > 0
                                            pp = polyfit(xTrC, yTrC, 1);
                                            if abs(pp(1)) > 1e-5
                                                pred = (yTe - pp(2)) / pp(1);
                                            else
                                                pred = yTe - (mean(yTrC) - mean(xTrC));
                                            end
                                        else
                                            pred = NaN;
                                        end

                                    case 3 % Proportional Correction
                                        vp = xTrC > 0 & yTrC > 0;
                                        if any(vp)
                                            ratio = mean(xTrC(vp) ./ yTrC(vp), 'omitnan');
                                            pred = yTe * ratio;
                                        else
                                            pred = NaN;
                                        end

                                    case 4 % Standard Deming (Fixed lambda = 1.0)
                                        [sl_d, ic_d] = fitWeightedDeming(app, xTrC, yTrC, 1.0);
                                        pred = (yTe - ic_d) / sl_d;

                                    case 5 % Linnet Weighted Deming (Tuned lambda strictly inside training fold)
                                        best_lam_fold = 1.0;
                                        best_rmse_wfold = inf;
                                        if numel(xTrC) >= 3
                                            for li = 1:numel(lam_grid)
                                                lv = lam_grid(li);
                                                [sl_w, ic_w] = fitWeightedDeming(app, xTrC, yTrC, lv);
                                                yh_w = (yTrC - ic_w) / sl_w;
                                                rmseW = sqrt(mean((yh_w - xTrC).^2, 'omitnan'));
                                                if rmseW < best_rmse_wfold
                                                    best_rmse_wfold = rmseW;
                                                    best_lam_fold = lv;
                                                end
                                            end
                                        end
                                        [sl_wd, ic_wd] = fitWeightedDeming(app, xTrC, yTrC, best_lam_fold);
                                        pred = (yTe - ic_wd) / sl_wd;

                                    case 6 % Passing-Bablok (Simplified)
                                        nT = numel(xTrC);
                                        slopesT = zeros(nT * (nT - 1) / 2, 1);
                                        idxT = 1;
                                        for ii = 1:nT-1
                                            for jj = ii+1:nT
                                                if xTrC(jj) ~= xTrC(ii)
                                                    slopesT(idxT) = (yTrC(jj) - yTrC(ii)) / (xTrC(jj) - xTrC(ii));
                                                else
                                                    slopesT(idxT) = NaN;
                                                end
                                                idxT = idxT + 1;
                                            end
                                        end
                                        slopesT = slopesT(~isnan(slopesT));
                                        if isempty(slopesT)
                                            pred = NaN;
                                        else
                                            slopeT = median(slopesT, 'omitnan');
                                            intT   = median(yTrC - slopeT .* xTrC, 'omitnan');
                                            if abs(slopeT) > 1e-5
                                                pred = (yTe - intT) / slopeT;
                                            else
                                                pred = yTe - (mean(yTrC) - mean(xTrC));
                                            end
                                        end

                                    case 7 % Hybrid Method (Tuned tau, W1, lambda strictly inside training fold)
                                        best_tau_fold = 0;
                                        best_w1_fold  = 8;
                                        best_lam_hfold = 1.0;
                                        best_rmse_hfold = inf;

                                        for wi = 1:nW1
                                            for ti = 1:nTau
                                                fp_tr = fastPairedFit(leaveIdx, ti, wi);
                                                cMaskH_tr = robustCleanMask(app, xTr, fp_tr);
                                                xTrH = xTr(cMaskH_tr);
                                                yTrH = fp_tr(cMaskH_tr);
                                                if numel(xTrH) < 2, continue; end

                                                for li = 1:numel(lam_grid)
                                                    lv = lam_grid(li);
                                                    [sl_h, ic_h] = fitWeightedDeming(app, xTrH, yTrH, lv);
                                                    yh_h = (yTrH - ic_h) / sl_h;
                                                    rmseH = sqrt(mean((yh_h - xTrH).^2, 'omitnan'));

                                                    corr_rough = fastRoughnessWin(ti, wi) / max(abs(sl_h), 1e-4);
                                                    if isnan(corr_rough), corr_rough = raw_rough; end
                                                    roughness_ratio = corr_rough / raw_rough;
                                                    pen = max(0, roughness_ratio - 1.0) * 0.25 * rmseH;

                                                    if rmseH + pen < best_rmse_hfold
                                                        best_rmse_hfold = rmseH + pen;
                                                        best_tau_fold   = tau_grid(ti);
                                                        best_w1_fold    = w1_grid(wi);
                                                        best_lam_hfold  = lv;
                                                    end
                                                end
                                            end
                                        end

                                        % Apply fold-selected hyperparameters to test sample i
                                        ti_fold = find(tau_grid == best_tau_fold, 1);
                                        wi_fold = find(w1_grid == best_w1_fold, 1);
                                        if isempty(ti_fold), ti_fold = 1; end
                                        if isempty(wi_fold), wi_fold = 2; end

                                        yFastTr = fastPairedFit(leaveIdx, ti_fold, wi_fold);
                                        yFastTe = fastPairedFit(i, ti_fold, wi_fold);
                                        cMaskH_f = robustCleanMask(app, xTr, yFastTr);
                                        [sl_hf, ic_hf] = fitWeightedDeming(app, xTr(cMaskH_f), yFastTr(cMaskH_f), best_lam_hfold);
                                        pred = (yFastTe - ic_hf) / sl_hf;
                                end
                                looErrors(i) = pred - xTe;
                            catch
                                looErrors(i) = NaN;
                            end
                        end
                        validErr = looErrors(~isnan(looErrors));
                        if numel(validErr) >= 2
                            candidateRMSE(ci)    = sqrt(mean(validErr.^2));
                            candidateLoASpan(ci) = 3.92 * std(validErr);
                        end
                    end

                    % Model Selection: Best LOO-CV RMSE with <1% LoA Span Tiebreaker
                    if n < 7
                        bestIdx = 1; 
                        bestName = candidateNames{1}; 
                        smallNWarning = true;
                    else
                        [sortedR, sortOrder] = sort(candidateRMSE);
                        if numel(sortOrder) >= 2 && sortedR(1) > 0 && ...
                                (sortedR(2) - sortedR(1)) / sortedR(1) < 0.01
                            tA2 = sortOrder(1); 
                            tB2 = sortOrder(2);
                            if candidateLoASpan(tB2) < candidateLoASpan(tA2)
                                bestIdx = tB2;
                            else
                                bestIdx = tA2;
                            end
                        else
                            bestIdx = sortOrder(1);
                        end
                        bestName = candidateNames{bestIdx}; 
                        smallNWarning = false;
                    end

                    % Fit deployment parameters for the winning model across the full fit-window dataset
                    cMaskFull = robustCleanMask(app, xABL_fit, yCDI_fit);
                    xCleanFull = xABL_fit(cMaskFull);
                    yCleanFull = yCDI_fit(cMaskFull);

                    switch bestIdx
                        case 1 % Bias
                            bias = mean(yCleanFull - xCleanFull, 'omitnan');
                            yCorrected = yCDI_fit - bias;
                            yCorrected(~cMaskFull) = NaN;
                            app.CorrectionModel.type = 'bias'; 
                            app.CorrectionModel.bias = bias;
                            app.CorrectionModel.formula = sprintf('Model: Raw_CDI = ABL %+.4f\nCorrected = Raw_CDI %+.4f\n(Fitted on %d/%d window pairs, %d robust [|diff - med| <= 4.5*MAD])', ...
                                bias, -bias, numel(xABL_fit), numel(xABL), sum(cMaskFull));

                        case 2 % OLS
                            pp = polyfit(xCleanFull, yCleanFull, 1);
                            if abs(pp(1)) < 1e-4, pp(1) = 1; pp(2) = mean(yCleanFull) - mean(xCleanFull); end
                            yCorrected = (yCDI_fit - pp(2)) / pp(1);
                            yCorrected(~cMaskFull) = NaN;
                            app.CorrectionModel.type = 'ols_abl_x'; 
                            app.CorrectionModel.slope = pp(1); 
                            app.CorrectionModel.intercept = pp(2);
                            app.CorrectionModel.formula = sprintf('Model: Raw_CDI = %.4f*ABL %+.4f\nCorrected = (Raw_CDI %+.4f) / %.4f\n(Fitted on %d/%d window pairs, %d robust [|diff - med| <= 4.5*MAD])', ...
                                pp(1), pp(2), -pp(2), pp(1), numel(xABL_fit), numel(xABL), sum(cMaskFull));

                        case 3 % Proportional
                            ratio = mean(xCleanFull ./ yCleanFull, 'omitnan');
                            yCorrected = yCDI_fit * ratio;
                            yCorrected(~cMaskFull) = NaN;
                            app.CorrectionModel.type = 'proportional'; 
                            app.CorrectionModel.ratio = ratio;
                            app.CorrectionModel.formula = sprintf('Ratio (ABL/CDI) = %.4f\nCorrected = Raw_CDI * %.4f\n(Fitted on %d/%d window pairs, %d robust [|diff - med| <= 4.5*MAD])', ...
                                ratio, ratio, numel(xABL_fit), numel(xABL), sum(cMaskFull));

                        case 4 % Deming (lambda=1.0)
                            [slope, intercept] = fitWeightedDeming(app, xCleanFull, yCleanFull, 1.0);
                            yCorrected = (yCDI_fit - intercept) / slope;
                            yCorrected(~cMaskFull) = NaN;
                            app.CorrectionModel.type = 'deming'; 
                            app.CorrectionModel.slope = slope; 
                            app.CorrectionModel.intercept = intercept; 
                            app.CorrectionModel.lam = 1.0;
                            app.DemingLambdaEditField.Value = 1.0;
                            app.CorrectionModel.formula = sprintf('Model: Raw_CDI = %.4f*ABL %+.4f  [Deming λ=1.00]\nCorrected = (Raw_CDI %+.4f) / %.4f\n(Fitted on %d/%d window pairs, %d robust [|diff - med| <= 4.5*MAD])', ...
                                slope, intercept, -intercept, slope, numel(xABL_fit), numel(xABL), sum(cMaskFull));

                        case 5 % Weighted Deming (Full tuning for deployment)
                            best_lam_deploy = 1.0;
                            best_rmse_wdeploy = inf;
                            for li = 1:numel(lam_grid)
                                lv = lam_grid(li);
                                [sl_w, ic_w] = fitWeightedDeming(app, xCleanFull, yCleanFull, lv);
                                yh_w = (yCleanFull - ic_w) / sl_w;
                                rmseW = sqrt(mean((yh_w - xCleanFull).^2, 'omitnan'));
                                if rmseW < best_rmse_wdeploy
                                    best_rmse_wdeploy = rmseW;
                                    best_lam_deploy = lv;
                                end
                            end
                            [slope, intercept] = fitWeightedDeming(app, xCleanFull, yCleanFull, best_lam_deploy);
                            yCorrected = (yCDI_fit - intercept) / slope;
                            yCorrected(~cMaskFull) = NaN;
                            app.CorrectionModel.type = 'weighted_deming'; 
                            app.CorrectionModel.slope = slope; 
                            app.CorrectionModel.intercept = intercept; 
                            app.CorrectionModel.lam = best_lam_deploy;
                            app.DemingLambdaEditField.Value = best_lam_deploy;
                            app.CorrectionModel.formula = sprintf('Model: Raw_CDI = %.4f*ABL %+.4f  [Linnet Weighted Deming λ=%.2f]\nCorrected = (Raw_CDI %+.4f) / %.4f\n(Fitted on %d/%d window pairs, %d robust [|diff - med| <= 4.5*MAD])', ...
                                slope, intercept, best_lam_deploy, -intercept, slope, numel(xABL_fit), numel(xABL), sum(cMaskFull));

                        case 6 % Passing-Bablok
                            n2 = numel(xCleanFull);
                            s2 = zeros(n2 * (n2 - 1) / 2, 1); 
                            i2 = 1;
                            for ii = 1:n2-1
                                for jj = ii+1:n2
                                    if xCleanFull(jj) ~= xCleanFull(ii)
                                        s2(i2) = (yCleanFull(jj) - yCleanFull(ii)) / (xCleanFull(jj) - xCleanFull(ii));
                                    else
                                        s2(i2) = NaN; 
                                    end
                                    i2 = i2 + 1;
                                end
                            end
                            s2 = s2(~isnan(s2));
                            slope = median(s2, 'omitnan');
                            if isnan(slope) || isinf(slope) || abs(slope) < 1e-4, slope = 1; end
                            intercept = median(yCleanFull - slope .* xCleanFull, 'omitnan');
                            if isnan(intercept), intercept = mean(yCleanFull) - mean(xCleanFull); end
                            yCorrected = (yCDI_fit - intercept) / slope;
                            yCorrected(~cMaskFull) = NaN;
                            app.CorrectionModel.type = 'passing-bablok'; 
                            app.CorrectionModel.slope = slope; 
                            app.CorrectionModel.intercept = intercept;
                            app.CorrectionModel.formula = sprintf('Model: Raw_CDI = %.4f*ABL %+.4f  [Simplified Passing-Bablok]\nCorrected = (Raw_CDI %+.4f) / %.4f\n(Fitted on %d/%d window pairs, %d robust [|diff - med| <= 4.5*MAD])', ...
                                slope, intercept, -intercept, slope, numel(xABL_fit), numel(xABL), sum(cMaskFull));

                        case 7 % Hybrid Method (Full tuning for deployment)
                            best_tau_dep = 0; best_w1_dep = 8; best_lam_dep = 1.0; best_rmse_dep = inf;
                            for wi = 1:nW1
                                for ti = 1:nTau
                                    fp_fit = fastPairedFit(:, ti, wi);
                                    cMaskH = robustCleanMask(app, xABL_fit, fp_fit);
                                    xCleanH = xABL_fit(cMaskH);
                                    yCleanH = fp_fit(cMaskH);
                                    if numel(xCleanH) < 2, continue; end
                                    
                                    for li = 1:numel(lam_grid)
                                        lv = lam_grid(li);
                                        [sl_h, ic_h] = fitWeightedDeming(app, xCleanH, yCleanH, lv);
                                        yh_h = (yCleanH - ic_h) / sl_h;
                                        rmseH = sqrt(mean((yh_h - xCleanH).^2, 'omitnan'));
                                        corr_rough = fastRoughnessWin(ti, wi) / max(abs(sl_h), 1e-4);
                                        if isnan(corr_rough), corr_rough = raw_rough; end
                                        roughness_ratio = corr_rough / raw_rough;
                                        pen = max(0, roughness_ratio - 1.0) * 0.25 * rmseH;
                                        if rmseH + pen < best_rmse_dep
                                            best_rmse_dep = rmseH + pen;
                                            best_tau_dep  = tau_grid(ti);
                                            best_w1_dep   = w1_grid(wi);
                                            best_lam_dep  = lv;
                                        end
                                    end
                                end
                            end

                            cdi_fast_fit = computeAsymmetricFastCDI(app, fullCDIVals, fullCDITime, best_w1_dep, best_tau_dep, best_tau_dep);
                            searchVals   = cdi_fast_fit(validCDIMask);
                            cdi_fast_paired = zeros(size(xABL));
                            for k = 1:numel(xABL)
                                [~, bi] = min(abs(searchTimes - pairedTimes(k)));
                                cdi_fast_paired(k) = searchVals(bi);
                            end
                            cMaskH = robustCleanMask(app, xABL_fit, cdi_fast_paired(fitWinMask));
                            xClean = xABL_fit(cMaskH); 
                            yClean = cdi_fast_paired(fitWinMask); 
                            yClean = yClean(cMaskH);

                            [slope, intercept] = fitWeightedDeming(app, xClean, yClean, best_lam_dep);
                            yCorrected = (cdi_fast_paired(fitWinMask) - intercept) / slope;
                            yCorrected(~cMaskH) = NaN;

                            app.CorrectionModel.type = 'hybrid';
                            app.CorrectionModel.tau_rise = best_tau_dep; 
                            app.CorrectionModel.tau_fall = best_tau_dep;
                            app.CorrectionModel.lam = best_lam_dep; 
                            app.CorrectionModel.slope = slope; 
                            app.CorrectionModel.intercept = intercept;
                            app.TauSpinner.Value = best_tau_dep; 
                            app.TauFallSpinner.Value = best_tau_dep;
                            app.SmoothW1Spinner.Value = best_w1_dep; 
                            app.DemingLambdaEditField.Value = best_lam_dep;
                            app.CorrectionModel.formula = sprintf('Step 1: CDI_s = movmean(CDI_raw, [%d 0])\nStep 2: dCDI/dt = movmean(diff(CDI_s)/dt, [%d 0])\nStep 3: CDI_fast = CDI_s + %.1f * (dCDI/dt)\nStep 4: CDI_fast = movmean(CDI_fast, [2 0])\nStep 5: CDI_corrected = (CDI_fast - %.4f) / %.4f  [Linnet λ=%.2f]\nHybrid filtering is causal: current output uses current and prior CDI samples only.\n(Fitted on %d/%d window pairs, %d robust [|diff - med| <= 4.5*MAD])', ...
                                max(0, best_w1_dep-1), max(0, best_w1_dep-1), best_tau_dep, intercept, slope, best_lam_dep, numel(xABL_fit), numel(xABL), sum(cMaskH));
                    end

                    [sortedRMSE, rankOrder] = sort(candidateRMSE);
                    rankText = cell(1, 7);
                    for ri = 1:7
                        marker = '';
                        if rankOrder(ri) == bestIdx, marker = '  ◀ winner'; end
                        rankText{ri} = sprintf('  #%d  %-35s  RMSE=%.4f  LoA span=%.4f%s', ...
                            ri, candidateNames{rankOrder(ri)}, sortedRMSE(ri), ...
                            candidateLoASpan(rankOrder(ri)), marker);
                    end

                    app.CorrectionModel.autoSelected  = true;
                    app.CorrectionModel.autoWinner    = bestName;
                    app.CorrectionModel.smallNWarning = smallNWarning;
                    app.CorrectionModel.autoRankText  = rankText;
                    app.CorrectionModel.autoRMSE      = candidateRMSE;
                    app.CorrectionModel.autoCandidates= candidateNames;

                otherwise
                    uialert(app.UIFigure, ['Unknown correction method selected: ', method], 'Method Error');
                    return;
            end
            
            if isempty(yCorrected) || all(isnan(yCorrected)) || all(isinf(yCorrected))
                uialert(app.UIFigure, 'Correction failed due to poor data (NaN/Inf produced).', 'Correction Error');
                return; 
            end
            
            yCorrected_display = yCorrected;
            app.CorrectionModel.yCorrected = yCorrected_display;
            app.CorrectionModel.timeVals = pairedTimes(fitWinMask);

            xABL_disp = xABL_fit;
            validFit = ~isnan(xABL_disp) & ~isnan(yCorrected_display) & ~isinf(yCorrected_display);
            if sum(validFit) >= 2
                app.CorrectionModel.pCorr = polyfit(xABL_disp(validFit), yCorrected_display(validFit), 1);
            else
                app.CorrectionModel.pCorr = [1, 0];
            end
            
            if std(xABL_disp,'omitnan')>0 && std(yCorrected_display,'omitnan')>0
                R_new = corrcoef(xABL_disp(validFit), yCorrected_display(validFit));
                app.CorrectionModel.r_new = R_new(1,2);
            else
                app.CorrectionModel.r_new = NaN;
            end

            diffNew  = yCorrected_display - xABL_disp;
            biasNew  = mean(diffNew, 'omitnan');
            sdNew    = std(diffNew,  'omitnan');
            
            % Compute BEFORE stats on the exact same valid (MAD-filtered) points
            diffBefore_clean = yCDI_fit(validFit) - xABL_disp(validFit);
            sdBefore_clean   = std(diffBefore_clean, 'omitnan');
            biasBefore_clean = mean(diffBefore_clean, 'omitnan');

            app.ImprovedBiasLabel.Text = sprintf('After Correction: Bias=%.4f', biasNew);
            app.ImprovedSDLabel.Text   = sprintf('SD=%.4f (on %d kept pairs)', sdNew, sum(validFit));

            loaBefore_clean = 3.92 * sdBefore_clean;
            loaAfter  = 3.92 * sdNew;
            sdReduction  = 100 * (1 - sdNew    / max(sdBefore_clean,  1e-10));
            loaReduction = 100 * (1 - loaAfter / max(loaBefore_clean, 1e-10));

            biasImproved = abs(biasNew) < abs(biasBefore_clean) - 1e-6;
            sdImproved   = sdReduction > 5;
            sdWorsened   = sdReduction < -5;

            if sdWorsened
                verdict = '✗ SD WORSENED';  vColor = [0.8 0.1 0.1];
            elseif biasImproved && sdImproved
                verdict = '✓ BIAS + SD IMPROVED';  vColor = [0.1 0.6 0.1];
            elseif biasImproved
                verdict = '~ BIAS ONLY';  vColor = [0.7 0.5 0.0];
            elseif sdImproved
                verdict = '~ SD ONLY';  vColor = [0.7 0.5 0.0];
            else
                verdict = '✗ NO IMPROVEMENT';  vColor = [0.8 0.1 0.1];
            end

            sdArrow  = '▼'; if sdReduction  < 0, sdArrow  = '▲'; end
            loaArrow = '▼'; if loaReduction < 0, loaArrow = '▲'; end
            app.CorrQualityLabel.Text = sprintf('%s  |  SD %s%.1f%%  LoA %s%.1f%%', ...
                verdict, sdArrow, abs(sdReduction), loaArrow, abs(loaReduction));
            app.CorrQualityLabel.FontColor = vColor;

            if isfield(app.CorrectionModel, 'formula')
                app.FormulaTextArea.Value = strsplit(app.CorrectionModel.formula, '\n');
            end

            app.CDI_Corrected_Table = app.CDI_Table;
            app.CDI_Corrected_Table.Time = app.CDI_Corrected_Table.Time + minutes(app.Stats.TimeShift);
            
            if ismember(param, app.CDI_Corrected_Table.Properties.VariableNames)
                origVals = app.CDI_Corrected_Table.(param);
                switch app.CorrectionModel.type
                    case 'hybrid'
                        cdi_fast = computeAsymmetricFastCDI(app, origVals, app.CDI_Corrected_Table.Time, ...
                            app.SmoothW1Spinner.Value, app.CorrectionModel.tau_rise, app.CorrectionModel.tau_fall);
                        if abs(app.CorrectionModel.slope) > 1e-10
                            full_corr = (cdi_fast - app.CorrectionModel.intercept) / app.CorrectionModel.slope;
                        else
                            full_corr = cdi_fast - app.CorrectionModel.intercept;
                        end
                        app.CDI_Corrected_Table.(param) = full_corr;
                    case 'bias'
                        app.CDI_Corrected_Table.(param) = origVals - app.CorrectionModel.bias;
                    case {'ols_abl_x', 'deming', 'weighted_deming', 'passing-bablok'}
                        app.CDI_Corrected_Table.(param) = (origVals - app.CorrectionModel.intercept) / app.CorrectionModel.slope;
                    case 'proportional'
                        app.CDI_Corrected_Table.(param) = origVals * app.CorrectionModel.ratio;
                end
                
                app.CorrectionModel.fullTime = app.CDI_Corrected_Table.Time;
                app.CorrectionModel.fullY = app.CDI_Corrected_Table.(param);
            end

            app.ShowCorrectedSwitch.Enable = 'on';
            app.ShowCorrectedSwitch.Value = true; 
            ShowCorrectedSwitchChanged(app, []);  
            
            app.KeepCorrectionSwitch.Enable = 'on';
            app.CorrViewSwitch.Enable = 'on';
            app.CorrViewSwitch.Value = true;
            CorrViewSwitchChanged(app, []); 
            
            app.BAViewSwitch.Enable = 'on';
            app.BAViewSwitch.Value = true;
            BAViewSwitchChanged(app, []); 
            
            if isfield(app.CorrectionModel, 'autoSelected') && app.CorrectionModel.autoSelected
                if isfield(app.CorrectionModel, 'smallNWarning') && app.CorrectionModel.smallNWarning
                    app.CorrectionStatusLabel.Text = sprintf('⚠ N=%d (<7): LOO-CV unreliable, defaulted to Bias Correction', numel(xABL_fit));
                    app.CorrectionStatusLabel.FontColor = [0.8 0 0];
                else
                    app.CorrectionStatusLabel.Text = sprintf('Auto → %s (RMSE=%.4f)', ...
                        app.CorrectionModel.autoWinner, min(app.CorrectionModel.autoRMSE));
                    app.CorrectionStatusLabel.FontColor = [0.2 0.6 0.2];
                end
            else
                app.CorrectionStatusLabel.Text = sprintf('%s applied', method);
                app.CorrectionStatusLabel.FontColor = [0.2 0.6 0.2];
            end
            
            try
                if app.FitWindowCheckBox.Value
                    fitStart = datetime(app.FitWindowStartEdit.Value, 'InputFormat', 'dd.MM.yyyy HH:mm');
                    fitEnd   = datetime(app.FitWindowEndEdit.Value,   'InputFormat', 'dd.MM.yyyy HH:mm');
                    [cvPct, nWin] = computeStabilityScore(app, fitStart, fitEnd, param, pairedTimes);
                    if isnan(cvPct)
                        app.StabilityScoreLabel.Text = sprintf('Window: %d pairs | CDI CV: N/A', nWin);
                    elseif cvPct < 5
                        app.StabilityScoreLabel.Text = sprintf('Window: %d pairs | CDI CV: %.1f%% ✓ stable', nWin, cvPct);
                        app.StabilityScoreLabel.FontColor = [0.1 0.6 0.1];
                    elseif cvPct < 15
                        app.StabilityScoreLabel.Text = sprintf('Window: %d pairs | CDI CV: %.1f%% ~ moderate', nWin, cvPct);
                        app.StabilityScoreLabel.FontColor = [0.7 0.5 0.0];
                    else
                        app.StabilityScoreLabel.Text = sprintf('Window: %d pairs | CDI CV: %.1f%% ⚠ unstable', nWin, cvPct);
                        app.StabilityScoreLabel.FontColor = [0.8 0.1 0.1];
                    end
                else
                    app.StabilityScoreLabel.Text = sprintf('Fitting on all %d pairs (no window)', numel(xABL));
                    app.StabilityScoreLabel.FontColor = [0.4 0.4 0.4];
                end
            catch
                app.StabilityScoreLabel.Text = '';
            end

            drawFitWindowShade(app);
            
            app.ExportCorrectedButton.Enable = 'on';
            app.ShowFormulaButton.Enable = 'on';
            app.ComparePlotsButton.Enable = 'on'; 
        end
        
        function ComparePlotsButtonPushed(app, ~)
            if isempty(app.CorrectionModel) || ~isfield(app.CorrectionModel, 'yCorrected')
                uialert(app.UIFigure, 'Please apply a correction first.', 'No Correction Data');
                return;
            end
            
            param = app.CurrentParam;
            xA = app.Stats.xABL;
            yOrig = app.Stats.yCDI;
            yCorr = app.CorrectionModel.yCorrected;
            
            % Ensure exact 1:1 pair comparison by removing MAD-filtered points from the baseline
            validPairs = ~isnan(yCorr) & ~isnan(xA);
            xA = xA(validPairs);
            yOrig = yOrig(validPairs);
            yCorr = yCorr(validPairs);
            
            screenSize = get(0, 'ScreenSize');
            winW = min(1200, screenSize(3) - 100);
            winH = min(900, screenSize(4) - 100);
            winX = max(50, round((screenSize(3) - winW) / 2));
            winY = max(50, round((screenSize(4) - winH) / 2));
            
            fig = uifigure('Name', sprintf('Before vs. After Analysis — %s', param), ...
                           'Position', [winX winY winW winH], ...
                           'Color', [0.95 0.95 0.95]);

            mFile = uimenu(fig, 'Text', 'File');
            uimenu(mFile, 'Text', 'Save Figure As...   (Ctrl+S)', 'Accelerator', 'S', ...
                'MenuSelectedFcn', @(~,~) ABL_CDI_Analyzer.saveCompareFigureAs(fig, param));
            uimenu(mFile, 'Text', 'Copy Figure to Clipboard', ...
                'MenuSelectedFcn', @(~,~) ABL_CDI_Analyzer.copyCompareFigure(fig));
            uimenu(mFile, 'Text', 'Print...', 'Separator', 'on', ...
                'MenuSelectedFcn', @(~,~) ABL_CDI_Analyzer.printCompareFigure(fig));
            uimenu(mFile, 'Text', 'Close Window', 'Separator', 'on', ...
                'MenuSelectedFcn', @(~,~) close(fig));

            mExport = uimenu(fig, 'Text', 'Export');
            uimenu(mExport, 'Text', 'Save as SVG (vector)...', ...
                'MenuSelectedFcn', @(~,~) ABL_CDI_Analyzer.exportCompareFigure(fig, param, 'svg'));
            uimenu(mExport, 'Text', 'Save as PNG (high resolution)...', ...
                'MenuSelectedFcn', @(~,~) ABL_CDI_Analyzer.exportCompareFigure(fig, param, 'png'));
            uimenu(mExport, 'Text', 'Save as PDF...', ...
                'MenuSelectedFcn', @(~,~) ABL_CDI_Analyzer.exportCompareFigure(fig, param, 'pdf'));
            uimenu(mExport, 'Text', 'Save as JPEG...', ...
                'MenuSelectedFcn', @(~,~) ABL_CDI_Analyzer.exportCompareFigure(fig, param, 'jpg'));
            uimenu(mExport, 'Text', 'Save as TIFF (300 DPI)...', 'Separator', 'on', ...
                'MenuSelectedFcn', @(~,~) ABL_CDI_Analyzer.exportCompareFigure(fig, param, 'tif'));

            gl = uigridlayout(fig, [2 2]);
            gl.RowHeight = {'1x', '1x'};
            gl.ColumnWidth = {'1x', '1x'};
            
            minValCorr = min([xA; yOrig; yCorr]); 
            maxValCorr = max([xA; yOrig; yCorr]);
            spanC = maxValCorr - minValCorr; if spanC==0, spanC=0.1; end
            padC = spanC * 0.1;
            limsCorr = [minValCorr-padC, maxValCorr+padC];
            
            diffOrig = yOrig - xA; 
            diffCorr = yCorr - xA;
            meanOrig = (yOrig + xA) / 2; 
            meanCorr = (yCorr + xA) / 2;
            
            minMean = min([meanOrig; meanCorr]); maxMean = max([meanOrig; meanCorr]);
            spanM = maxMean - minMean; if spanM==0, spanM=0.1; end
            padM = spanM * 0.1; limsMean = [minMean-padM, maxMean+padM];
            
            minDiff = min([diffOrig; diffCorr]); maxDiff = max([diffOrig; diffCorr]);
            spanD = maxDiff - minDiff; if spanD==0, spanD=0.1; end
            padD = spanD * 0.1; limsDiff = [minDiff-padD, maxDiff+padD];
            
            ax1 = uiaxes(gl); ax1.Layout.Row = 1; ax1.Layout.Column = 1; ax1.Box = 'on';
            hold(ax1, 'on');
            scatter(ax1, xA, yOrig, 50, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);
            plot(ax1, limsCorr, limsCorr, 'k--', 'DisplayName', 'Identity');
            if ~any(isnan(app.Stats.p))
                xF = linspace(limsCorr(1), limsCorr(2), 100);
                plot(ax1, xF, polyval(app.Stats.p, xF), 'r-', 'LineWidth', 1.5);
            end
            title(ax1, sprintf('ORIGINAL Correlation (r=%.3f)', app.Stats.r));
            xlabel(ax1, ['ABL ' param]); ylabel(ax1, ['CDI ' param]);
            xlim(ax1, limsCorr); ylim(ax1, limsCorr); grid(ax1, 'on'); hold(ax1, 'off');
            
            ax2 = uiaxes(gl); ax2.Layout.Row = 1; ax2.Layout.Column = 2; ax2.Box = 'on';
            hold(ax2, 'on');
            scatter(ax2, xA, yCorr, 50, 'filled', 'MarkerFaceColor', [0.2 0.7 0.2]);
            plot(ax2, limsCorr, limsCorr, 'k--', 'DisplayName', 'Identity');
            if ~any(isnan(app.CorrectionModel.pCorr))
                xF = linspace(limsCorr(1), limsCorr(2), 100);
                plot(ax2, xF, polyval(app.CorrectionModel.pCorr, xF), 'r-', 'LineWidth', 1.5);
            end
            title(ax2, sprintf('CORRECTED Correlation (r=%.3f)', app.CorrectionModel.r_new));
            xlabel(ax2, ['ABL ' param]); ylabel(ax2, ['Corrected CDI ' param]);
            xlim(ax2, limsCorr); ylim(ax2, limsCorr); grid(ax2, 'on'); hold(ax2, 'off');
            
            ax3 = uiaxes(gl); ax3.Layout.Row = 2; ax3.Layout.Column = 1; ax3.Box = 'on';
            hold(ax3, 'on');
            biasO = mean(diffOrig, 'omitnan'); sdO = std(diffOrig, 'omitnan');
            scatter(ax3, meanOrig, diffOrig, 50, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);
            plot(ax3, limsMean, [biasO biasO], 'b-', 'LineWidth', 1.5);
            plot(ax3, limsMean, [biasO+1.96*sdO biasO+1.96*sdO], 'r--', 'LineWidth', 1);
            plot(ax3, limsMean, [biasO-1.96*sdO biasO-1.96*sdO], 'r--', 'LineWidth', 1);
            title(ax3, sprintf('ORIGINAL Bland-Altman (Bias: %.2f)', biasO));
            xlabel(ax3, 'Mean (ABL, CDI)'); ylabel(ax3, 'Diff (CDI - ABL)');
            xlim(ax3, limsMean); ylim(ax3, limsDiff); grid(ax3, 'on'); hold(ax3, 'off');
            
            ax4 = uiaxes(gl); ax4.Layout.Row = 2; ax4.Layout.Column = 2; ax4.Box = 'on';
            hold(ax4, 'on');
            biasC = mean(diffCorr, 'omitnan'); sdC = std(diffCorr, 'omitnan');
            scatter(ax4, meanCorr, diffCorr, 50, 'filled', 'MarkerFaceColor', [0.2 0.7 0.2]);
            plot(ax4, limsMean, [biasC biasC], 'b-', 'LineWidth', 1.5);
            plot(ax4, limsMean, [biasC+1.96*sdC biasC+1.96*sdC], 'g--', 'LineWidth', 1);
            plot(ax4, limsMean, [biasC-1.96*sdC biasC-1.96*sdC], 'g--', 'LineWidth', 1);
            title(ax4, sprintf('CORRECTED Bland-Altman (Bias: %.2f)', biasC));
            xlabel(ax4, 'Mean (ABL, CDI_c_o_r_r)'); ylabel(ax4, 'Diff (CDI_c_o_r_r - ABL)');
            xlim(ax4, limsMean); ylim(ax4, limsDiff); grid(ax4, 'on'); hold(ax4, 'off');
        end

        function ShowFormulaButtonPushed(app, ~)
            if isempty(app.CorrectionModel) || ~isfield(app.CorrectionModel, 'formula')
                return;
            end
            
            param = app.CurrentParam;
            mdl   = app.CorrectionModel;
            isAuto = isfield(mdl, 'autoSelected') && mdl.autoSelected;
            
            xA = app.Stats.xABL;
            yC = app.Stats.yCDI;
            
            latexStr = {};
            latexStr{end+1} = '% ----- GENERAL THEORETICAL EQUATION -----';
            
            switch mdl.type
                case 'hybrid'
                    yCorr = mdl.yCorrected;
                    w_sm = max(0, app.SmoothW1Spinner.Value - 1);
                    latexStr{end+1} = sprintf('$$ \\text{CDI}_{s} = \\text{movmean}(\\text{CDI}_{\\text{raw}},\\, %d) $$', w_sm);
                    latexStr{end+1} = sprintf('$$ \\frac{d\\text{CDI}}{dt} = \\text{movmean}\\!\\left(\\frac{\\Delta \\text{CDI}_s}{\\Delta t},\\, %d\\right) $$', w_sm);
                    latexStr{end+1} = '$$ \text{CDI}_{\text{fast}} = \text{CDI}_{s} + \tau_{\text{dir}} \cdot \frac{d\text{CDI}}{dt} $$';
                    latexStr{end+1} = '$$ \text{CDI}_{\text{fast}} = \text{movmean}(\text{CDI}_{\text{fast}},\, 2) $$';
                    latexStr{end+1} = '$$ \text{CDI}_{\text{corrected}} = \frac{\text{CDI}_{\text{fast}} - \beta_0}{\beta_1} \quad [\text{Weighted Deming}] $$';
                    latexStr{end+1} = '';
                    latexStr{end+1} = '% ----- SPECIFIC APPLICATION EQUATION -----';
                    latexStr{end+1} = sprintf('$$ \\tau_{\\text{rise}} = %.1f \\text{ min}, \\quad \\tau_{\\text{fall}} = %.1f \\text{ min}, \\quad \\beta_0 = %.4f, \\quad \\beta_1 = %.4f, \\quad \\lambda = %.2f $$', mdl.tau_rise, mdl.tau_fall, mdl.intercept, mdl.slope, mdl.lam);
                case 'weighted_deming'
                    yCorr = (yC - mdl.intercept) / mdl.slope;
                    latexStr{end+1} = '$$ \text{CDI}_{\text{corrected}} = \frac{\text{CDI}_{\text{raw}} - \beta_0}{\beta_1} \quad [\text{Linnet Weighted Deming: } w_i = 1/\hat{u}_i^2] $$';
                    latexStr{end+1} = '';
                    latexStr{end+1} = '% ----- SPECIFIC APPLICATION EQUATION -----';
                    latexStr{end+1} = sprintf('$$ \\text{CDI}_{\\text{corrected}} = \\frac{\\text{CDI}_{\\text{raw}} %+.4f}{%.4f}, \\quad \\lambda = %.2f $$', -mdl.intercept, mdl.slope, mdl.lam);
                case 'bias'
                    yCorr = yC - mdl.bias;
                    latexStr{end+1} = '$$ \text{CDI}_{\text{corrected}} = \text{CDI}_{\text{raw}} - \text{Bias} $$';
                    latexStr{end+1} = '';
                    latexStr{end+1} = '% ----- SPECIFIC APPLICATION EQUATION -----';
                    if mdl.bias > 0
                        latexStr{end+1} = sprintf('$$ \\text{CDI}_{\\text{corrected}} = \\text{CDI}_{\\text{raw}} - %.4f $$', mdl.bias);
                    else
                        latexStr{end+1} = sprintf('$$ \\text{CDI}_{\\text{corrected}} = \\text{CDI}_{\\text{raw}} + %.4f $$', abs(mdl.bias));
                    end
                case 'ols_abl_x'
                    yCorr = (yC - mdl.intercept) / mdl.slope;
                    latexStr{end+1} = '$$ \text{CDI}_{\text{corrected}} = \frac{\text{CDI}_{\text{raw}} - \beta_0}{\beta_1} $$';
                    latexStr{end+1} = '';
                    latexStr{end+1} = '% ----- SPECIFIC APPLICATION EQUATION -----';
                    latexStr{end+1} = sprintf('$$ \\text{CDI}_{\\text{corrected}} = \\frac{\\text{CDI}_{\\text{raw}} %+.4f}{%.4f} $$', -mdl.intercept, mdl.slope);
                case 'proportional'
                    yCorr = yC * mdl.ratio;
                    latexStr{end+1} = '$$ \text{CDI}_{\text{corrected}} = \text{CDI}_{\text{raw}} \times \text{Ratio} $$';
                    latexStr{end+1} = '';
                    latexStr{end+1} = '% ----- SPECIFIC APPLICATION EQUATION -----';
                    latexStr{end+1} = sprintf('$$ \\text{CDI}_{\\text{corrected}} = \\text{CDI}_{\\text{raw}} \\times %.4f $$', mdl.ratio);
                case 'deming'
                    yCorr = (yC - mdl.intercept) / mdl.slope;
                    latexStr{end+1} = '$$ \text{CDI}_{\text{corrected}} = \frac{\text{CDI}_{\text{raw}} - \beta_0}{\beta_1} $$';
                    latexStr{end+1} = '';
                    latexStr{end+1} = '% ----- SPECIFIC APPLICATION EQUATION -----';
                    latexStr{end+1} = sprintf('$$ \\text{CDI}_{\\text{corrected}} = \\frac{\\text{CDI}_{\\text{raw}} %+.4f}{%.4f}, \\quad \\lambda = %.2f $$', -mdl.intercept, mdl.slope, mdl.lam);
                case 'passing-bablok'
                    yCorr = (yC - mdl.intercept) / mdl.slope;
                    latexStr{end+1} = '$$ \text{CDI}_{\text{corrected}} = \frac{\text{CDI}_{\text{raw}} - \beta_0}{\beta_1} $$';
                    latexStr{end+1} = '';
                    latexStr{end+1} = '% ----- SPECIFIC APPLICATION EQUATION -----';
                    latexStr{end+1} = sprintf('$$ \\text{CDI}_{\\text{corrected}} = \\frac{\\text{CDI}_{\\text{raw}} %+.4f}{%.4f} $$', -mdl.intercept, mdl.slope);
            end

            % Mask arrays to compare only pairs surviving the MAD filter
            validPairs = ~isnan(yCorr) & ~isnan(xA);
            dBefore = yC(validPairs) - xA(validPairs);
            dAfter  = yCorr(validPairs) - xA(validPairs);
            nBeforeMAD = numel(xA);
            nAfterMAD  = sum(validPairs);
            
            mName = mdl.type; mName(1) = upper(mName(1));
            mName = strrep(mName, '_', ' ');
            if isAuto, mName = [mName ' (auto-selected)']; end
            
            howTo = '';
            switch mdl.type
                case 'hybrid'
                    w1_val = app.SmoothW1Spinner.Value;
                    w_sm = max(0, w1_val - 1);
                    howTo = sprintf('THE CONCEPT:\nDynamic response filtering + Linnet Weighted Deming regression.\nHybrid filtering is causal: current output uses current and prior CDI samples only.\n\nSTEP 1 - PRE-SMOOTHING (W1-1=%d):\nCDI_s = movmean(CDI_raw, [%d 0])\n\nSTEP 2 - DERIVATIVE LEAD:\nCDI_fast = CDI_s + tau * (dCDI_s/dt)\n\nSTEP 3 - LINNET WEIGHTED DEMING:\nFitted with lambda=%.2f using 1/u^2 variance weighting.\nCDI_corrected = (CDI_fast - %.4f) / %.4f', w_sm, w_sm, mdl.lam, mdl.intercept, mdl.slope);
                case 'weighted_deming'
                    howTo = sprintf('THE ALGEBRA:\nModel: Raw_CDI = (%.4f * ABL) %+.4f\nCorrected = (Raw_CDI %+.4f) / %.4f\n\nHOW IT WAS CALCULATED:\n- Linnet Weighted Deming Regression.\n- Weights w_i = 1 / u_hat_i^2 account for proportional error variance.\n- Lambda = %.2f variance ratio.', mdl.slope, mdl.intercept, -mdl.intercept, mdl.slope, mdl.lam);
                case 'bias'
                    if mdl.bias > 0
                        howTo = sprintf('THE ALGEBRA:\nModel: Raw_CDI = ABL + %.4f\nCorrected = Raw_CDI - %.4f\n\nHOW IT WAS CALCULATED:\n- Constant offset shift based on robust median difference.', mdl.bias, mdl.bias);
                    else
                        howTo = sprintf('THE ALGEBRA:\nModel: Raw_CDI = ABL - %.4f\nCorrected = Raw_CDI + %.4f\n\nHOW IT WAS CALCULATED:\n- Constant offset shift based on robust median difference.', abs(mdl.bias), abs(mdl.bias));
                    end
                case 'ols_abl_x'
                    howTo = sprintf('THE ALGEBRA:\nModel: Raw_CDI = (%.4f * ABL) %+.4f\nCorrected = (Raw_CDI %+.4f) / %.4f', mdl.slope, mdl.intercept, -mdl.intercept, mdl.slope);
                case 'proportional'
                    howTo = sprintf('THE ALGEBRA:\nRatio = Mean(ABL / Raw_CDI) = %.4f\nCorrected = Raw_CDI * %.4f', mdl.ratio, mdl.ratio);
                case 'deming'
                    howTo = sprintf('THE ALGEBRA:\nModel: Raw_CDI = (%.4f * ABL) %+.4f  [Deming λ=%.2f]\nCorrected = (Raw_CDI %+.4f) / %.4f', mdl.slope, mdl.intercept, mdl.lam, -mdl.intercept, mdl.slope);
                case 'passing-bablok'
                    howTo = sprintf('THE ALGEBRA:\nModel: Raw_CDI = (%.4f * ABL) %+.4f  [Passing-Bablok]\nCorrected = (Raw_CDI %+.4f) / %.4f', mdl.slope, mdl.intercept, -mdl.intercept, mdl.slope);
            end
            
            screenSize = get(0, 'ScreenSize');
            winW = min(950, screenSize(3) - 100);
            if isAuto, winH = min(900, screenSize(4) - 50);
            else, winH = min(800, screenSize(4) - 50); end
            winX = max(50, round((screenSize(3) - winW) / 2));
            winY = max(50, round((screenSize(4) - winH) / 2));
            
            fw = uifigure('Name', sprintf('Correction Report — %s', param), ...
                          'Position', [winX winY winW winH], ...
                          'Color', [1 1 1], 'Resize', 'on');
            
            if isAuto
                numRows = 5;
                rowHeights = {28, '3x', 'fit', '2x', '1x'};
            else
                numRows = 4;
                rowHeights = {28, '3x', 'fit', '1x'};
            end
            
            mg = uigridlayout(fw, [numRows 2]);
            mg.RowHeight = rowHeights;
            mg.ColumnWidth = {'1x', '1x'};
            mg.Padding = [10 10 10 10];
            mg.RowSpacing = 6;
            mg.ColumnSpacing = 8;
            
            titleLbl = uilabel(mg);
            titleLbl.Layout.Row = 1; titleLbl.Layout.Column = [1 2];
            titleLbl.Text = sprintf('  %s  |  %s  |  N = %d', param, mName, app.Stats.N);
            titleLbl.FontSize = 14; titleLbl.FontWeight = 'bold';
            titleLbl.FontColor = [1 1 1];
            titleLbl.BackgroundColor = [0.2 0.35 0.6];
            titleLbl.HorizontalAlignment = 'left';
            
            leftPanel = uipanel(mg, 'Title', 'Formula & Application', ...
                'FontWeight', 'bold', 'FontSize', 11, ...
                'BackgroundColor', [0.97 1.0 0.97]);
            leftPanel.Layout.Row = 2; leftPanel.Layout.Column = 1;
            
            lg = uigridlayout(leftPanel, [2 1]);
            lg.RowHeight = {'fit', '1x'};
            lg.Padding = [6 4 6 4]; lg.RowSpacing = 4;
            
            fBox = uitextarea(lg);
            fBox.Layout.Row = 1;
            fBox.FontSize = 14; 
            fBox.FontName = 'Consolas';
            fBox.FontWeight = 'bold';
            fBox.Editable = 'off';
            fBox.BackgroundColor = [0.93 0.98 0.93];
            fBox.Value = {mdl.formula};
            
            hBox = uitextarea(lg);
            hBox.Layout.Row = 2;
            hBox.FontSize = 12; 
            hBox.FontName = 'Consolas';
            hBox.Editable = 'off';
            hBox.BackgroundColor = [0.97 0.97 1.0];
            hBox.Value = strsplit(howTo, '\n');
            
            rightPanel = uipanel(mg, 'Title', 'Before vs After Correction', ...
                'FontWeight', 'bold', 'FontSize', 11, ...
                'BackgroundColor', [1.0 0.98 0.94]);
            rightPanel.Layout.Row = 2; rightPanel.Layout.Column = 2;
            
            rg = uigridlayout(rightPanel, [2 1]);
            rg.RowHeight = {'1x', 'fit'};
            rg.Padding = [6 4 6 4];
            
            sBox = uitextarea(rg);
            sBox.Layout.Row = 1;
            sBox.FontSize = 11; 
            sBox.FontName = 'Consolas';
            sBox.Editable = 'off';
            sBox.BackgroundColor = [1.0 0.98 0.94];
            
            statsTextLines = {
                sprintf('                BEFORE         AFTER')
                sprintf('Bias:        %+9.4f     %+9.4f', mean(dBefore, 'omitnan'), mean(dAfter, 'omitnan'))
                sprintf('SD:           %9.4f      %9.4f', std(dBefore, 'omitnan'), std(dAfter, 'omitnan'))
                sprintf('95%% LoA:   [%+.2f, %+.2f]  [%+.2f, %+.2f]', ...
                    mean(dBefore,'omitnan')-1.96*std(dBefore,'omitnan'), mean(dBefore,'omitnan')+1.96*std(dBefore,'omitnan'), ...
                    mean(dAfter,'omitnan')-1.96*std(dAfter,'omitnan'), mean(dAfter,'omitnan')+1.96*std(dAfter,'omitnan'))
                sprintf('Max |err|:    %9.4f      %9.4f', max(abs(dBefore)), max(abs(dAfter)))
                ''
                sprintf('Bias reduction:  %.1f%%', 100*(1 - abs(mean(dAfter,'omitnan'))/max(abs(mean(dBefore,'omitnan')),1e-10)))
                sprintf('SD change:       %+.1f%%', 100*(std(dAfter,'omitnan')/max(std(dBefore,'omitnan'),1e-10) - 1))
                ''
                sprintf('MAD Filter: |diff - median| > 4.5*MAD')
                sprintf('Excluded from stats. Kept %d / %d pairs.', nAfterMAD, nBeforeMAD)
            };
            sBox.Value = statsTextLines;
            
            expStatsBtn = uibutton(rg, 'push');
            expStatsBtn.Layout.Row = 2;
            expStatsBtn.Text = '💾 Export as SVG';
            expStatsBtn.FontSize = 9;
            expStatsBtn.FontWeight = 'bold';
            expStatsBtn.BackgroundColor = [0.15 0.45 0.75];
            expStatsBtn.FontColor = [1 1 1];
            expStatsBtn.ButtonPushedFcn = @(~,~) ABL_CDI_Analyzer.exportStatsAsSVG(statsTextLines, param, mName);
            
            interpPanel = uipanel(mg, 'Title', 'Interpretation', ...
                'FontWeight', 'bold', 'FontSize', 11, ...
                'BackgroundColor', [0.96 0.96 0.96]);
            interpPanel.Layout.Row = 3; interpPanel.Layout.Column = [1 2];
            
            ig = uigridlayout(interpPanel, [1 1]);
            ig.Padding = [6 3 6 3];
            
            interpLines = {};
            biasAfter = abs(mean(dAfter, 'omitnan'));
            sdAfter = std(dAfter, 'omitnan');
            if biasAfter < 0.5 && sdAfter < 2
                interpLines{end+1} = 'Bias + SD improved - residual bias near zero, tight agreement.';
            elseif biasAfter < 1.0
                interpLines{end+1} = 'Bias reduced, scatter remains.';
            else
                interpLines{end+1} = 'Limited improvement - evaluate non-linearities or extreme sensor drift.';
            end
            
            iBox = uilabel(ig);
            iBox.Layout.Row = 1;
            iBox.Text = strjoin(interpLines, '  ');
            iBox.FontSize = 11; 
            iBox.WordWrap = 'on';
            
            if isAuto
                autoPanel = uipanel(mg, 'Title', 'Auto-Selection Report (7-Model LOO Cross-Validation)', ...
                    'FontWeight', 'bold', 'FontSize', 11, ...
                    'BackgroundColor', [0.93 0.96 1.0]);
                autoPanel.Layout.Row = 4; autoPanel.Layout.Column = [1 2];
                
                ag = uigridlayout(autoPanel, [1 1]);
                ag.Padding = [6 4 6 4];
                
                autoLines = {};
                autoLines{end+1} = '   Primary criterion: Fold-wise LOO-CV RMSE (parameters tuned strictly inside folds)';
                autoLines{end+1} = '   Final parameters refitted on the complete selected fit-window after ranking.';
                autoLines{end+1} = '';
                for ri = 1:numel(mdl.autoRankText)
                    autoLines{end+1} = mdl.autoRankText{ri}; 
                end
                autoLines{end+1} = '';
                autoLines{end+1} = sprintf('Winner: %s', mdl.autoWinner);
                
                aBox = uitextarea(ag);
                aBox.Layout.Row = 1;
                aBox.FontSize = 11; 
                aBox.FontName = 'Consolas';
                aBox.Editable = 'off';
                aBox.BackgroundColor = [0.93 0.96 1.0];
                aBox.Value = autoLines;
            end
            
            latexPanel = uipanel(mg, 'Title', 'LaTeX Code for Equations', ...
                'FontWeight', 'bold', 'FontSize', 11, ...
                'BackgroundColor', [0.95 0.95 0.95]);
            latexPanel.Layout.Row = numRows; latexPanel.Layout.Column = [1 2];
            
            latexGrid = uigridlayout(latexPanel, [1 1]);
            latexGrid.Padding = [6 4 6 4];
            
            lBox = uitextarea(latexGrid);
            lBox.Layout.Row = 1;
            lBox.FontSize = 11;
            lBox.FontName = 'Consolas';
            lBox.Editable = 'on'; 
            lBox.Value = latexStr;
        end
        
        function ExportCorrectedButtonPushed(app, ~)
            if isempty(app.CDI_Corrected_Table)
                uialert(app.UIFigure, 'No corrected data available. Apply correction first.', 'No Data');
                return;
            end
            [file, path] = uiputfile({'*.csv';'*.xlsx'}, 'Save Corrected CDI Data As');
            if isequal(file,0), return; end
            fullpath = fullfile(path, file);
            try
                writetable(app.CDI_Corrected_Table, fullpath);
                app.StatusLabel.Text = ['Corrected CDI exported to: ' file];
            catch ME
                uialert(app.UIFigure, ['Export failed: ' ME.message], 'Export Error');
            end
        end

        function ShowCorrectedSwitchChanged(app, ~)
            if ~isempty(app.CorrectedTrendLine) && isvalid(app.CorrectedTrendLine)
                delete(app.CorrectedTrendLine);
                app.CorrectedTrendLine = [];
            end

            if app.ShowCorrectedSwitch.Value && isfield(app.CorrectionModel, 'fullY') && isfield(app.CorrectionModel, 'fullTime')
                hold(app.TimeAxes, 'on');
                t = app.CorrectionModel.fullTime;
                y = app.CorrectionModel.fullY;

                corrLabel = app.CorrectionMethodDropDown.Value;
                if isfield(app.CorrectionModel, 'autoSelected') && app.CorrectionModel.autoSelected
                    corrLabel = ['Auto: ' app.CorrectionModel.autoWinner];
                end
                corrLabel = ['CDI corrected (' corrLabel ')'];
                if app.FitWindowCheckBox.Value
                    corrLabel = [corrLabel ' [window coeffs applied to full recording]'];
                end

                app.CorrectedTrendLine = plot(app.TimeAxes, t, y, ...
                    '-', 'DisplayName', corrLabel, ...
                    'LineWidth', 1.6, 'Color', [0.15 0.65 0.15]);
                
                legend(app.TimeAxes, 'Location', 'best');
                hold(app.TimeAxes, 'off');
            end
        end

        function ShowOriginalCDISwitchChanged(app, ~)
            if ~isempty(app.OriginalCDILine) && isvalid(app.OriginalCDILine)
                if app.ShowOriginalCDISwitch.Value
                    app.OriginalCDILine.Visible = 'on';
                else
                    app.OriginalCDILine.Visible = 'off';
                end
                legend(app.TimeAxes, 'Location', 'best');
            end
        end

        function KeepCorrectionSwitchChanged(app, ~)
            if app.KeepCorrectionSwitch.Value
                if isfield(app.CorrectionModel, 'fullY') && isfield(app.CorrectionModel, 'fullTime')
                    methodLabel = app.CorrectionMethodDropDown.Value;
                    if isfield(app.CorrectionModel, 'autoSelected') && app.CorrectionModel.autoSelected
                        methodLabel = ['Auto: ' app.CorrectionModel.autoWinner];
                    end
                    
                    if ~isempty(app.KeptCorrectionLines)
                        removeIdx = [];
                        for ki = 1:numel(app.KeptCorrectionLines)
                            if isvalid(app.KeptCorrectionLines{ki}) && ...
                                    strcmp(app.KeptCorrectionLines{ki}.DisplayName, methodLabel)
                                delete(app.KeptCorrectionLines{ki});
                                removeIdx(end+1) = ki; 
                            end
                        end
                        app.KeptCorrectionLines(removeIdx) = [];
                        app.KeptCorrectionData(removeIdx) = [];
                    end
                    
                    app.KeptCorrectionData{end+1} = struct( ...
                        'name', methodLabel, ...
                        'timeVals', app.CorrectionModel.fullTime, ...
                        'yVals', app.CorrectionModel.fullY);
                    
                    hold(app.TimeAxes, 'on');
                    keptColors = [0.8 0.6 0; 0.5 0 0.5; 0 0.5 0.5; 0.8 0.2 0.2; 0.2 0.2 0.8];
                    ci = mod(numel(app.KeptCorrectionData)-1, size(keptColors,1)) + 1;
                    kLine = plot(app.TimeAxes, ...
                        app.CorrectionModel.fullTime, app.CorrectionModel.fullY, ...
                        '--', 'DisplayName', methodLabel, ...
                        'LineWidth', 1.3, 'Color', keptColors(ci,:));
                    if isempty(app.KeptCorrectionLines)
                        app.KeptCorrectionLines = {kLine};
                    else
                        app.KeptCorrectionLines{end+1} = kLine;
                    end
                    legend(app.TimeAxes, 'Location', 'best');
                    hold(app.TimeAxes, 'off');
                end
            else
                if ~isempty(app.KeptCorrectionLines)
                    for ki = 1:numel(app.KeptCorrectionLines)
                        if isvalid(app.KeptCorrectionLines{ki})
                            delete(app.KeptCorrectionLines{ki});
                        end
                    end
                    app.KeptCorrectionLines = {};
                end
                app.KeptCorrectionData = {};
                legend(app.TimeAxes, 'Location', 'best');
            end
        end

        % --- UI INITIALIZATION ---
        function createComponents(app)
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1440 850];
            app.UIFigure.Name = 'ABL-CDI Analyzer';
            app.UIFigure.Color = [0.95 0.95 0.95];

            app.GridLayout = uigridlayout(app.UIFigure, [1 2]);
            app.GridLayout.ColumnWidth = {310, '1x'};
            app.GridLayout.RowHeight = {'1x'};

            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;
            app.LeftPanel.Title = 'Controls & Calibration';
            app.LeftPanel.FontWeight = 'bold';
            app.LeftPanel.Scrollable = 'on';
            app.LeftPanel.AutoResizeChildren = 'off';  

            app.LeftGrid = uigridlayout(app.LeftPanel, [22 1]);
            app.LeftGrid.ColumnWidth = {'1x'};
            app.LeftGrid.RowHeight = repmat({'fit'}, 1, 22);
            app.LeftGrid.Padding = [10 10 10 10];
            app.LeftGrid.RowSpacing = 6;
            app.LeftGrid.Scrollable = 'on';

            app.LoadABLButton = uibutton(app.LeftGrid, 'push');
            app.LoadABLButton.Layout.Row = 1;
            app.LoadABLButton.Text = '1️⃣ Load ABL Reference Data';
            app.LoadABLButton.FontSize = 11;
            app.LoadABLButton.ButtonPushedFcn = createCallbackFcn(app, @LoadABLButtonPushed, true);

            app.LoadCDIButton = uibutton(app.LeftGrid, 'push');
            app.LoadCDIButton.Layout.Row = 2;
            app.LoadCDIButton.Text = '2️⃣ Load CDI Continuous Log';
            app.LoadCDIButton.FontSize = 11;
            app.LoadCDIButton.ButtonPushedFcn = createCallbackFcn(app, @LoadCDIButtonPushed, true);

            app.PatientIDLabel = uilabel(app.LeftGrid);
            app.PatientIDLabel.Layout.Row = 3;
            app.PatientIDLabel.Text = 'Patient ID:';
            app.PatientIDLabel.FontWeight = 'bold';
            app.PatientIDLabel.FontSize = 10;

            app.PatientIDDropDown = uidropdown(app.LeftGrid);
            app.PatientIDDropDown.Layout.Row = 4;
            app.PatientIDDropDown.Items = {'All Patients'};
            app.PatientIDDropDown.Value = 'All Patients';
            app.PatientIDDropDown.FontSize = 11;
            app.PatientIDDropDown.ValueChangedFcn = createCallbackFcn(app, @PatientIDDropDownValueChanged, true);

            app.ParamLabel = uilabel(app.LeftGrid);
            app.ParamLabel.Layout.Row = 5;
            app.ParamLabel.Text = 'Select Analyte Parameter:';
            app.ParamLabel.FontWeight = 'bold';
            app.ParamLabel.FontSize = 10;

            app.ParamDropDown = uidropdown(app.LeftGrid);
            app.ParamDropDown.Layout.Row = 6;
            app.ParamDropDown.Items = {'Load files first...'};
            app.ParamDropDown.FontSize = 11;

            app.TimeToleranceLabel = uilabel(app.LeftGrid);
            app.TimeToleranceLabel.Layout.Row = 7;
            app.TimeToleranceLabel.Text = 'Matching Tolerance (min):';
            app.TimeToleranceLabel.FontWeight = 'bold';
            app.TimeToleranceLabel.FontSize = 10;

            app.TimeToleranceSpinner = uispinner(app.LeftGrid);
            app.TimeToleranceSpinner.Layout.Row = 8;
            app.TimeToleranceSpinner.Limits = [0.5 60];
            app.TimeToleranceSpinner.Value = 5;
            app.TimeToleranceSpinner.Step = 0.5; 
            app.TimeToleranceSpinner.ValueChangedFcn = createCallbackFcn(app, @TimeToleranceSpinnerValueChanged, true);

            app.TimeShiftLabel = uilabel(app.LeftGrid);
            app.TimeShiftLabel.Layout.Row = 9;
            app.TimeShiftLabel.Text = 'CDI Timeline Shift (min):';
            app.TimeShiftLabel.FontWeight = 'bold';
            app.TimeShiftLabel.FontSize = 10;

            app.TimeShiftSpinner = uispinner(app.LeftGrid);
            app.TimeShiftSpinner.Layout.Row = 10;
            app.TimeShiftSpinner.Limits = [-1440 1440]; 
            app.TimeShiftSpinner.Value = 0; 
            app.TimeShiftSpinner.Step = 1;
            app.TimeShiftSpinner.ValueChangedFcn = createCallbackFcn(app, @TimeShiftSpinnerValueChanged, true);
            
            app.AutoShiftButton = uibutton(app.LeftGrid, 'push');
            app.AutoShiftButton.Layout.Row = 11;
            app.AutoShiftButton.Text = '⏱️ Auto-Detect Timeline Shift';
            app.AutoShiftButton.FontSize = 10;
            app.AutoShiftButton.FontWeight = 'bold';
            app.AutoShiftButton.BackgroundColor = [0.9 0.9 0.9];
            app.AutoShiftButton.ButtonPushedFcn = createCallbackFcn(app, @AutoShiftButtonPushed, true);

            app.AnalyzeButton = uibutton(app.LeftGrid, 'push');
            app.AnalyzeButton.Layout.Row = 12;
            app.AnalyzeButton.Text = '3️⃣ RUN ANALYSIS';
            app.AnalyzeButton.FontSize = 12;
            app.AnalyzeButton.FontWeight = 'bold';
            app.AnalyzeButton.BackgroundColor = [0.3 0.75 0.93];
            app.AnalyzeButton.FontColor = [1 1 1];
            app.AnalyzeButton.ButtonPushedFcn = createCallbackFcn(app, @AnalyzeButtonPushed, true);

            % Fitting Window Panel
            app.FitWindowPanel = uipanel(app.LeftGrid);
            app.FitWindowPanel.Layout.Row = 13;
            app.FitWindowPanel.Title = '📐 Fitting Window';
            app.FitWindowPanel.FontWeight = 'bold';
            app.FitWindowPanel.FontSize = 10;
            app.FitWindowPanel.BackgroundColor = [0.94 0.97 1.0];
            app.FitWindowPanel.Scrollable = 'off';

            fwGrid = uigridlayout(app.FitWindowPanel, [8 1]);
            fwGrid.ColumnWidth = {'1x'};
            fwGrid.RowHeight = {'fit','fit','fit','fit','fit','fit','fit','fit'};
            fwGrid.Padding = [8 6 8 6];
            fwGrid.RowSpacing = 4;

            app.FitWindowCheckBox = uicheckbox(fwGrid);
            app.FitWindowCheckBox.Layout.Row = 1;
            app.FitWindowCheckBox.Text = 'Enable fitting window';
            app.FitWindowCheckBox.Value = false;
            app.FitWindowCheckBox.FontSize = 9;
            app.FitWindowCheckBox.FontWeight = 'bold';
            app.FitWindowCheckBox.ValueChangedFcn = createCallbackFcn(app, @FitWindowCheckBoxChanged, true);

            app.FitWindowStartLabel = uilabel(fwGrid);
            app.FitWindowStartLabel.Layout.Row = 2;
            app.FitWindowStartLabel.Text = 'Window start (dd.MM.yyyy HH:mm):';
            app.FitWindowStartLabel.FontSize = 8;

            app.FitWindowStartEdit = uieditfield(fwGrid, 'text');
            app.FitWindowStartEdit.Layout.Row = 3;
            app.FitWindowStartEdit.Value = '';
            app.FitWindowStartEdit.FontSize = 9;
            app.FitWindowStartEdit.FontName = 'Consolas';
            app.FitWindowStartEdit.Placeholder = 'e.g. 01.02.2026 11:00';

            app.FitWindowEndLabel = uilabel(fwGrid);
            app.FitWindowEndLabel.Layout.Row = 4;
            app.FitWindowEndLabel.Text = 'Window end (dd.MM.yyyy HH:mm):';
            app.FitWindowEndLabel.FontSize = 8;

            app.FitWindowEndEdit = uieditfield(fwGrid, 'text');
            app.FitWindowEndEdit.Layout.Row = 5;
            app.FitWindowEndEdit.Value = '';
            app.FitWindowEndEdit.FontSize = 9;
            app.FitWindowEndEdit.FontName = 'Consolas';
            app.FitWindowEndEdit.Placeholder = 'e.g. 01.02.2026 18:00';

            app.FitWindowAutoButton = uibutton(fwGrid, 'push');
            app.FitWindowAutoButton.Layout.Row = 6;
            app.FitWindowAutoButton.Text = '🎯 Auto Window Overlap';
            app.FitWindowAutoButton.FontSize = 9;
            app.FitWindowAutoButton.FontWeight = 'bold';
            app.FitWindowAutoButton.BackgroundColor = [0.15 0.55 0.35];
            app.FitWindowAutoButton.FontColor = [1 1 1];
            app.FitWindowAutoButton.ButtonPushedFcn = createCallbackFcn(app, @FitWindowAutoButtonPushed, true);

            applyBtn = uibutton(fwGrid, 'push');
            applyBtn.Layout.Row = 7;
            applyBtn.Text = '✓ Apply typed window';
            applyBtn.FontSize = 9;
            applyBtn.FontWeight = 'bold';
            applyBtn.BackgroundColor = [0.25 0.55 0.85];
            applyBtn.FontColor = [1 1 1];
            applyBtn.ButtonPushedFcn = createCallbackFcn(app, @FitWindowApplyButtonPushed, true);

            app.StabilityScoreLabel = uilabel(fwGrid);
            app.StabilityScoreLabel.Layout.Row = 8;
            app.StabilityScoreLabel.Text = 'Apply correction to see stability';
            app.StabilityScoreLabel.FontSize = 8;
            app.StabilityScoreLabel.FontColor = [0.4 0.4 0.4];
            app.StabilityScoreLabel.WordWrap = 'on';  

            % Stats Panel
            app.StatsPanel = uipanel(app.LeftGrid);
            app.StatsPanel.Layout.Row = 14;
            app.StatsPanel.Title = 'Descriptive Agreement Statistics';
            app.StatsPanel.FontWeight = 'bold';
            app.StatsPanel.FontSize = 10;
            app.StatsPanel.BackgroundColor = [1 1 1];

            app.StatsGrid = uigridlayout(app.StatsPanel, [10 1]);
            app.StatsGrid.ColumnWidth = {'1x'};
            app.StatsGrid.RowHeight = repmat({15}, 1, 10);
            app.StatsGrid.Padding = [5 2 5 2];
            app.StatsGrid.RowSpacing = 1;

            app.NPairsLabel = uilabel(app.StatsGrid);
            app.NPairsLabel.Layout.Row = 1;
            app.NPairsLabel.Text = 'N Pairs: --';
            app.NPairsLabel.FontSize = 9;

            app.BiasLabel = uilabel(app.StatsGrid);
            app.BiasLabel.Layout.Row = 2;
            app.BiasLabel.Text = 'Bias: --';
            app.BiasLabel.FontSize = 9;

            app.SDLabel = uilabel(app.StatsGrid);
            app.SDLabel.Layout.Row = 3;
            app.SDLabel.Text = 'SD: --';
            app.SDLabel.FontSize = 9;

            app.LOALabel = uilabel(app.StatsGrid);
            app.LOALabel.Layout.Row = 4;
            app.LOALabel.Text = '95% LoA: --';
            app.LOALabel.FontSize = 9;

            app.CorrelationLabel = uilabel(app.StatsGrid);
            app.CorrelationLabel.Layout.Row = 5;
            app.CorrelationLabel.Text = 'r = --';
            app.CorrelationLabel.FontSize = 9;

            app.RegressionLabel = uilabel(app.StatsGrid);
            app.RegressionLabel.Layout.Row = 6;
            app.RegressionLabel.Text = 'Regression: --';
            app.RegressionLabel.FontSize = 9;

            app.ImprovedBiasLabel = uilabel(app.StatsGrid);
            app.ImprovedBiasLabel.Layout.Row = 7;
            app.ImprovedBiasLabel.Text = 'After Correction: --';
            app.ImprovedBiasLabel.FontSize = 9;
            app.ImprovedBiasLabel.FontWeight = 'bold';

            app.ImprovedSDLabel = uilabel(app.StatsGrid);
            app.ImprovedSDLabel.Layout.Row = 8;
            app.ImprovedSDLabel.Text = '';
            app.ImprovedSDLabel.FontSize = 9;
            app.ImprovedSDLabel.FontWeight = 'bold';

            app.CorrQualityLabel = uilabel(app.StatsGrid);
            app.CorrQualityLabel.Layout.Row = 9;
            app.CorrQualityLabel.Text = '';
            app.CorrQualityLabel.FontSize = 9;
            app.CorrQualityLabel.FontWeight = 'bold';
            app.CorrQualityLabel.FontColor = [0.2 0.6 0.2];

            app.SlopeWarningLabel = uilabel(app.StatsGrid);
            app.SlopeWarningLabel.Layout.Row = 10;
            app.SlopeWarningLabel.Text = '';
            app.SlopeWarningLabel.FontSize = 8;
            app.SlopeWarningLabel.FontColor = [0.8 0.1 0.1];
            app.SlopeWarningLabel.WordWrap = 'on';
            app.SlopeWarningLabel.Visible = 'off';

            % Correction Panel
            app.CorrectionPanel = uipanel(app.LeftGrid);
            app.CorrectionPanel.Layout.Row = 15;
            app.CorrectionPanel.Title = '🔧 Method Correction & Tuning';
            app.CorrectionPanel.FontWeight = 'bold';
            app.CorrectionPanel.FontSize = 10;
            app.CorrectionPanel.BackgroundColor = [1 1 0.95];
            app.CorrectionPanel.Scrollable = 'on';
            app.CorrectionPanel.AutoResizeChildren = 'off';

            app.CorrectionGrid = uigridlayout(app.CorrectionPanel, [19 1]); 
            app.CorrectionGrid.ColumnWidth = {'1x'};
            app.CorrectionGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 140, 'fit'};
            app.CorrectionGrid.Padding = [8 8 8 8];
            app.CorrectionGrid.RowSpacing = 4;
            app.CorrectionGrid.Scrollable = 'on';

            app.CorrectionMethodLabel = uilabel(app.CorrectionGrid);
            app.CorrectionMethodLabel.Layout.Row = 1;
            app.CorrectionMethodLabel.Text = 'Correction Model:';
            app.CorrectionMethodLabel.FontSize = 9;
            app.CorrectionMethodLabel.FontWeight = 'bold';

            app.CorrectionMethodDropDown = uidropdown(app.CorrectionGrid);
            app.CorrectionMethodDropDown.Layout.Row = 2;
            app.CorrectionMethodDropDown.Items = {'Auto (Best Model)', 'Weighted Deming (Linnet)', 'Hybrid (Time-Series + Deming)', 'Deming Regression', 'Bias Correction', 'OLS (ABL is X)', 'Proportional Correction', 'Passing-Bablok'};
            app.CorrectionMethodDropDown.Value = 'Auto (Best Model)';
            app.CorrectionMethodDropDown.FontSize = 10;

            app.DemingLambdaLabel = uilabel(app.CorrectionGrid);
            app.DemingLambdaLabel.Layout.Row = 3;
            app.DemingLambdaLabel.Text = 'Deming λ (Variance Ratio ABL/CDI):';
            app.DemingLambdaLabel.FontSize = 9;

            app.DemingLambdaEditField = uieditfield(app.CorrectionGrid, 'numeric');
            app.DemingLambdaEditField.Layout.Row = 4;
            app.DemingLambdaEditField.Value = 1.0;
            app.DemingLambdaEditField.Limits = [0.001 1000];

            app.TauLabel = uilabel(app.CorrectionGrid);
            app.TauLabel.Layout.Row = 5;
            app.TauLabel.Text = 'Hybrid τ_rise (min):';
            app.TauLabel.FontSize = 9;

            app.TauSpinner = uispinner(app.CorrectionGrid);
            app.TauSpinner.Layout.Row = 6;
            app.TauSpinner.Value = 2.0;
            app.TauSpinner.Step = 0.5;
            app.TauSpinner.Limits = [0 20];

            app.TauFallLabel = uilabel(app.CorrectionGrid);
            app.TauFallLabel.Layout.Row = 7;
            app.TauFallLabel.Text = 'Hybrid τ_fall (min):';
            app.TauFallLabel.FontSize = 9;

            app.TauFallSpinner = uispinner(app.CorrectionGrid);
            app.TauFallSpinner.Layout.Row = 8;
            app.TauFallSpinner.Value = 1.5;
            app.TauFallSpinner.Step = 0.5;
            app.TauFallSpinner.Limits = [0 20];

            app.SmoothW1Label = uilabel(app.CorrectionGrid);
            app.SmoothW1Label.Layout.Row = 9;
            app.SmoothW1Label.Text = 'Smooth window W1 (samples):';
            app.SmoothW1Label.FontSize = 9;

            app.SmoothW1Spinner = uispinner(app.CorrectionGrid);
            app.SmoothW1Spinner.Layout.Row = 10;
            app.SmoothW1Spinner.Value = 8;
            app.SmoothW1Spinner.Step = 4;
            app.SmoothW1Spinner.Limits = [2 64];

            app.AutoTuneButton = uibutton(app.CorrectionGrid, 'push');
            app.AutoTuneButton.Layout.Row = 11;
            app.AutoTuneButton.Text = '⚙️ Auto-Tune Params (Single Fit)';
            app.AutoTuneButton.FontSize = 9;
            app.AutoTuneButton.FontWeight = 'bold';
            app.AutoTuneButton.BackgroundColor = [0.9 0.9 0.9];
            app.AutoTuneButton.ButtonPushedFcn = createCallbackFcn(app, @AutoTuneButtonPushed, true);

            app.ApplyCorrectionButton = uibutton(app.CorrectionGrid, 'push');
            app.ApplyCorrectionButton.Layout.Row = 12;
            app.ApplyCorrectionButton.Text = '✓ Apply Correction';
            app.ApplyCorrectionButton.FontSize = 11;
            app.ApplyCorrectionButton.FontWeight = 'bold';
            app.ApplyCorrectionButton.BackgroundColor = [0.4 0.8 0.4];
            app.ApplyCorrectionButton.FontColor = [1 1 1];
            app.ApplyCorrectionButton.Enable = 'off';
            app.ApplyCorrectionButton.ButtonPushedFcn = createCallbackFcn(app, @ApplyCorrectionButtonPushed, true);

            app.ExportCorrectedButton = uibutton(app.CorrectionGrid, 'push');
            app.ExportCorrectedButton.Layout.Row = 13;
            app.ExportCorrectedButton.Text = '💾 Export Corrected';
            app.ExportCorrectedButton.FontSize = 10;
            app.ExportCorrectedButton.FontWeight = 'bold';
            app.ExportCorrectedButton.BackgroundColor = [0.15 0.45 0.75];
            app.ExportCorrectedButton.FontColor = [1 1 1];
            app.ExportCorrectedButton.Enable = 'off';
            app.ExportCorrectedButton.ButtonPushedFcn = createCallbackFcn(app, @ExportCorrectedButtonPushed, true);

            app.ShowFormulaButton = uibutton(app.CorrectionGrid, 'push');
            app.ShowFormulaButton.Layout.Row = 14;
            app.ShowFormulaButton.Text = '📐 Show Formula & LaTeX';
            app.ShowFormulaButton.FontSize = 10;
            app.ShowFormulaButton.FontWeight = 'bold';
            app.ShowFormulaButton.BackgroundColor = [0.95 0.7 0.3];
            app.ShowFormulaButton.Enable = 'off';
            app.ShowFormulaButton.ButtonPushedFcn = createCallbackFcn(app, @ShowFormulaButtonPushed, true);

            app.ComparePlotsButton = uibutton(app.CorrectionGrid, 'push');
            app.ComparePlotsButton.Layout.Row = 15;
            app.ComparePlotsButton.Text = '📊 Compare Before/After';
            app.ComparePlotsButton.FontSize = 10;
            app.ComparePlotsButton.FontWeight = 'bold';
            app.ComparePlotsButton.BackgroundColor = [0.6 0.4 0.8];
            app.ComparePlotsButton.FontColor = [1 1 1];
            app.ComparePlotsButton.Enable = 'off';
            app.ComparePlotsButton.ButtonPushedFcn = createCallbackFcn(app, @ComparePlotsButtonPushed, true);

            app.CorrectionStatusLabel = uilabel(app.CorrectionGrid);
            app.CorrectionStatusLabel.Layout.Row = 16;
            app.CorrectionStatusLabel.Text = 'No correction applied';
            app.CorrectionStatusLabel.FontSize = 8;
            app.CorrectionStatusLabel.FontColor = [0.5 0.5 0.5];
            app.CorrectionStatusLabel.HorizontalAlignment = 'center';

            app.FormulaLabel = uilabel(app.CorrectionGrid);
            app.FormulaLabel.Layout.Row = 17;
            app.FormulaLabel.Text = 'Mathematical Model Formula:';
            app.FormulaLabel.FontSize = 9;
            app.FormulaLabel.FontWeight = 'bold';

            app.FormulaTextArea = uitextarea(app.CorrectionGrid);
            app.FormulaTextArea.Layout.Row = 18;
            app.FormulaTextArea.FontSize = 9;
            app.FormulaTextArea.FontName = 'Consolas';
            app.FormulaTextArea.Editable = 'off';
            app.FormulaTextArea.Value = {'Apply correction, then click "Show Formula"'};
            app.FormulaTextArea.BackgroundColor = [0.98 0.98 1];
            app.FormulaTextArea.WordWrap = 'on';

            app.SmallNWarningLabel = uilabel(app.CorrectionGrid);
            app.SmallNWarningLabel.Layout.Row = 19;
            app.SmallNWarningLabel.Text = '';
            app.SmallNWarningLabel.FontSize = 9;
            app.SmallNWarningLabel.FontWeight = 'bold';
            app.SmallNWarningLabel.FontColor = [0.8 0 0];
            app.SmallNWarningLabel.HorizontalAlignment = 'left';
            app.SmallNWarningLabel.WordWrap = 'on';
            app.SmallNWarningLabel.Visible = 'off';

            app.ExportButton = uibutton(app.LeftGrid);
            app.ExportButton.Layout.Row = 16;
            app.ExportButton.Text = '📊 Export Analysis Stats';
            app.ExportButton.FontSize = 10;
            app.ExportButton.Enable = 'off';
            app.ExportButton.ButtonPushedFcn = createCallbackFcn(app, @ExportButtonPushed, true);

            app.ExportFigureButton = uibutton(app.LeftGrid);
            app.ExportFigureButton.Layout.Row = 17;
            app.ExportFigureButton.Text = '🖼️ Export Report Figures (SVG)';
            app.ExportFigureButton.FontSize = 10;
            app.ExportFigureButton.Enable = 'off';
            app.ExportFigureButton.BackgroundColor = [0.25 0.55 0.35];
            app.ExportFigureButton.FontColor = [1 1 1];
            app.ExportFigureButton.ButtonPushedFcn = createCallbackFcn(app, @ExportFigureButtonPushed, true);

            app.StatusLabel = uilabel(app.LeftGrid);
            app.StatusLabel.Layout.Row = 18;
            app.StatusLabel.Text = 'Ready';
            app.StatusLabel.FontColor = [0.2 0.6 0.2];
            app.StatusLabel.FontSize = 10;

            % Right Visualization Panel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;
            app.RightPanel.Title = 'Diagnostic Visualization & Analysis';
            app.RightPanel.FontWeight = 'bold';

            app.RightGrid = uigridlayout(app.RightPanel, [4 3]);
            app.RightGrid.ColumnWidth = {'1x', '1x', '1x'};
            app.RightGrid.RowHeight = {'1.2x', 'fit', 'fit', '1x'};
            app.RightGrid.Padding = [10 10 10 10];

            app.TimeAxes = uiaxes(app.RightGrid);
            app.TimeAxes.Layout.Row = 1;
            app.TimeAxes.Layout.Column = [1 3];
            app.TimeAxes.Box = 'on';

            app.ShowOriginalCDISwitch = uicheckbox(app.RightGrid);
            app.ShowOriginalCDISwitch.Layout.Row = 2;
            app.ShowOriginalCDISwitch.Layout.Column = 1;
            app.ShowOriginalCDISwitch.Text = 'Show original CDI';
            app.ShowOriginalCDISwitch.Value = true;
            app.ShowOriginalCDISwitch.FontSize = 10;
            app.ShowOriginalCDISwitch.ValueChangedFcn = createCallbackFcn(app, @ShowOriginalCDISwitchChanged, true);

            app.ShowCorrectedSwitch = uicheckbox(app.RightGrid);
            app.ShowCorrectedSwitch.Layout.Row = 2;
            app.ShowCorrectedSwitch.Layout.Column = 2;
            app.ShowCorrectedSwitch.Text = 'Show corrected CDI';
            app.ShowCorrectedSwitch.FontSize = 10;
            app.ShowCorrectedSwitch.Enable = 'off';
            app.ShowCorrectedSwitch.ValueChangedFcn = createCallbackFcn(app, @ShowCorrectedSwitchChanged, true);

            app.KeepCorrectionSwitch = uicheckbox(app.RightGrid);
            app.KeepCorrectionSwitch.Layout.Row = 2;
            app.KeepCorrectionSwitch.Layout.Column = 3;
            app.KeepCorrectionSwitch.Text = 'Pin current correction';
            app.KeepCorrectionSwitch.FontSize = 10;
            app.KeepCorrectionSwitch.Enable = 'off';
            app.KeepCorrectionSwitch.ValueChangedFcn = createCallbackFcn(app, @KeepCorrectionSwitchChanged, true);

            app.CorrViewSwitch = uicheckbox(app.RightGrid);
            app.CorrViewSwitch.Layout.Row = 3;
            app.CorrViewSwitch.Layout.Column = 1;
            app.CorrViewSwitch.Text = 'Show corrected scatter';
            app.CorrViewSwitch.FontSize = 10;
            app.CorrViewSwitch.Enable = 'off';
            app.CorrViewSwitch.ValueChangedFcn = createCallbackFcn(app, @CorrViewSwitchChanged, true);

            app.BAViewSwitch = uicheckbox(app.RightGrid);
            app.BAViewSwitch.Layout.Row = 3;
            app.BAViewSwitch.Layout.Column = 2;
            app.BAViewSwitch.Text = 'Show corrected Bland-Altman';
            app.BAViewSwitch.FontSize = 10;
            app.BAViewSwitch.Enable = 'off';
            app.BAViewSwitch.ValueChangedFcn = createCallbackFcn(app, @BAViewSwitchChanged, true);

            app.CorrelationAxes = uiaxes(app.RightGrid);
            app.CorrelationAxes.Layout.Row = 4;
            app.CorrelationAxes.Layout.Column = 1;
            app.CorrelationAxes.Box = 'on';

            app.BlandAltmanAxes = uiaxes(app.RightGrid);
            app.BlandAltmanAxes.Layout.Row = 4;
            app.BlandAltmanAxes.Layout.Column = [2 3];
            app.BlandAltmanAxes.Box = 'on';

            app.UIFigure.Visible = 'on';
        end
    end

    methods (Static)
        function exportStatsAsSVG(statsLines, param, methodName)
            [file, path] = uiputfile({'*.svg','SVG vector (*.svg)'}, ...
                'Save Statistics As', sprintf('%s_%s_statistics.svg', param, strrep(methodName,' ','_')));
            if isequal(file,0), return; end
            fullpath = fullfile(path, file);

            hF = figure('Visible','off','Units','centimeters','Position',[0 0 10 8]);
            set(hF,'Color','white');
            ax = axes(hF,'Visible','off');
            ax.Units='normalized'; ax.Position=[0 0 1 1];

            nLines = numel(statsLines);
            yStep = 0.88 / (nLines + 2);

            text(ax, 0.05, 0.94, sprintf('Before vs After - %s  (%s)', param, methodName), ...
                'Units','normalized','FontSize',9,'FontWeight','bold',...
                'Color',[0.15 0.25 0.55],'Interpreter','none','FontName','Consolas');

            for si = 1:nLines
                yPos = 0.88 - (si-1)*yStep;
                txt = statsLines{si};
                if isempty(txt), continue; end
                isSummary = contains(txt,'reduction') || contains(txt,'change');
                fw = 'normal'; col = [0.1 0.1 0.1];
                if isSummary, col = [0.15 0.45 0.15]; fw = 'bold'; end
                text(ax, 0.05, yPos, txt, ...
                    'Units','normalized','FontSize',9,'FontName','Consolas',...
                    'FontWeight',fw,'Color',col,'Interpreter','none');
            end

            rectangle('Parent',ax,'Position',[0.01 0.01 0.98 0.98],...
                'EdgeColor',[0.5 0.5 0.5],'LineWidth',0.8);

            print(hF, fullpath, '-dsvg');
            close(hF);
        end

        function exportCompareFigure(fig, param, fmt)
            safeParam = regexprep(string(param), '[^a-zA-Z0-9_]', '_');
            defaultName = sprintf('Compare_Before_After_%s', safeParam);
            
            switch lower(fmt)
                case 'svg'
                    spec = {'*.svg', 'SVG vector image (*.svg)'}; ext = '.svg';
                case 'png'
                    spec = {'*.png', 'PNG image, 300 DPI (*.png)'}; ext = '.png';
                case 'pdf'
                    spec = {'*.pdf', 'PDF document (*.pdf)'}; ext = '.pdf';
                case {'jpg','jpeg'}
                    spec = {'*.jpg', 'JPEG image, 300 DPI (*.jpg)'}; ext = '.jpg';
                case {'tif','tiff'}
                    spec = {'*.tif', 'TIFF image, 300 DPI (*.tif)'}; ext = '.tif';
                otherwise
                    spec = {'*.png', 'PNG image (*.png)'}; ext = '.png';
            end
            
            [file, path] = uiputfile(spec, 'Save Compare Figure As', defaultName);
            if isequal(file, 0), return; end
            fullpath = fullfile(path, file);
            
            try
                if strcmpi(ext, '.svg') || strcmpi(ext, '.pdf')
                    exportgraphics(fig, fullpath, 'ContentType', 'vector');
                elseif strcmpi(ext, '.tif')
                    exportgraphics(fig, fullpath, 'Resolution', 300);
                else
                    exportgraphics(fig, fullpath, 'Resolution', 300);
                end
                msgbox(sprintf('Saved to:\n%s', fullpath), 'Export complete', 'help');
            catch ME
                uialert(fig, ['Export failed: ' ME.message], 'Export Error');
            end
        end
        
        function saveCompareFigureAs(fig, param)
            choice = uiconfirm(fig, 'Choose output format:', 'Save Figure As', ...
                'Options', {'SVG (vector)', 'PNG (300 DPI)', 'PDF', 'Cancel'}, ...
                'DefaultOption', 1, 'CancelOption', 4);
            switch choice
                case 'SVG (vector)'
                    ABL_CDI_Analyzer.exportCompareFigure(fig, param, 'svg');
                case 'PNG (300 DPI)'
                    ABL_CDI_Analyzer.exportCompareFigure(fig, param, 'png');
                case 'PDF'
                    ABL_CDI_Analyzer.exportCompareFigure(fig, param, 'pdf');
            end
        end
        
        function copyCompareFigure(fig)
            try
                copygraphics(fig, 'Resolution', 300, 'BackgroundColor', 'white');
                disp('Compare figure copied to clipboard.');
            catch ME
                uialert(fig, ['Copy to clipboard failed: ' ME.message], 'Copy Error');
            end
        end
        
        function printCompareFigure(fig)
            try
                print(fig, '-dprinter');
            catch ME
                uialert(fig, ['Print failed: ' ME.message], 'Print Error');
            end
        end
    end

    methods (Access = public)
        function app = ABL_CDI_Analyzer
            createComponents(app)
            registerApp(app, app.UIFigure)
        end

        function delete(app)
            if isvalid(app.UIFigure)
                delete(app.UIFigure)
            end
        end
    end
end