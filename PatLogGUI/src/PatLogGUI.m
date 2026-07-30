function PatLogGUI()

    % Helper function to convert any value to numeric
    function numVal = toNumber(val)
        if isnumeric(val)
            numVal = val;
        elseif ischar(val) || isstring(val)
            valStr = char(val);
            valStr = strrep(valStr, ',', '.');
            numVal = str2double(valStr);
        else
            numVal = NaN;
        end
    end
    
    % Helper function to get parameter name with units
    function paramWithUnit = getParamWithUnit(paramName)
        switch paramName
            case 'pH'
                paramWithUnit = 'pH';
            case 'pO2'
                paramWithUnit = 'pO2 (mmHg)';
            case 'pCO2'
                paramWithUnit = 'pCO2 (mmHg)';
            case 'T (°C)'
                paramWithUnit = 'T (°C)';
            case 'tHb'
                paramWithUnit = 'tHb (g/dL)';
            case 'FIO2'
                paramWithUnit = 'FIO2 (%)';
            case 'sO2'
                paramWithUnit = 'sO2 (%)';
            case 'K+'
                paramWithUnit = 'K+ (mmol/L)';
            case 'Lac'
                paramWithUnit = 'Lac (mmol/L)';
            case 'Na+'
                paramWithUnit = 'Na+ (mmol/L)';
            case 'Cl-'
                paramWithUnit = 'Cl- (mmol/L)';
            otherwise
                paramWithUnit = paramName;
        end
    end
    
    % Helper function to get selected patient IDs (supports multi-selection)
    function [patientNums, patientLabel] = getSelectedPatients()
        multiIDText = strtrim(multiPatientField.Value);
        
        if ~isempty(multiIDText)
            % Multi-patient mode
            idStrings = strsplit(multiIDText, ',');
            patientNums = [];
            for i = 1:length(idStrings)
                idStr = strtrim(idStrings{i});
                num = str2double(idStr);
                if ~isnan(num)
                    patientNums(end+1) = num;
                end
            end
            if ~isempty(patientNums)
                patientLabel = strjoin(arrayfun(@num2str, patientNums, 'UniformOutput', false), ',');
            else
                patientNums = [];
                patientLabel = '';
            end
        else
            % Single patient mode from dropdown
            selectedPatient = patientDropdown.Value;
            if strcmp(selectedPatient, 'All Patients')
                patientNums = [];  % Empty means all
                patientLabel = 'All Patients';
            else
                patientNum = str2double(selectedPatient);
                if ~isnan(patientNum)
                    patientNums = patientNum;
                    patientLabel = selectedPatient;
                else
                    patientNums = [];
                    patientLabel = '';
                end
            end
        end
    end

    %% Initialize GUI variables
    rawData = [];
    cleanData = [];
    headers = [];
    paramNames = {'pH', 'pO2', 'pCO2', 'T (°C)', 'tHb', 'FIO2', 'sO2', 'K+', 'Lac', 'Na+', 'Cl-'};
    paramIndices = [];
    patientIDCol = 0;
    currentFile = '';
    
    %% Create main figure
    fig = uifigure('Name', 'PatLog Blood Gas Analyzer - Clinical Data Management', ...
        'Position', [50 50 1500 850]);
    
    % Create main grid layout
    mainGrid = uigridlayout(fig, [1 2]);
    mainGrid.ColumnWidth = {300, '1x'};
    
    %% LEFT PANEL - Controls
    leftPanel = uipanel(mainGrid, 'Title', 'Clinical Controls', ...
        'FontWeight', 'bold', 'FontSize', 11, ...
        'Scrollable', 'on', ...
        'AutoResizeChildren', 'off');
    leftGrid = uigridlayout(leftPanel, [27 1]);
    leftGrid.RowHeight = repmat({'fit'}, 1, 27);
    leftGrid.Padding = [15 15 15 15];
    leftGrid.RowSpacing = 10;
    leftGrid.Scrollable = 'on';
    
    uilabel(leftGrid, 'Text', 'Blood Gas Analysis System', ...
        'FontSize', 18, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'FontColor', [0.2 0.4 0.7]);
    
    uilabel(leftGrid, 'Text', '──── Data Import ────', ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        'FontColor', [0.3 0.3 0.3]);
    
    selectFileBtn = uibutton(leftGrid, 'Text', '📁 Import Patient Data File', ...
        'ButtonPushedFcn', @(btn,event) selectFile(), ...
        'BackgroundColor', [0.25 0.55 0.85], 'FontColor', 'white', ...
        'FontSize', 11, 'FontWeight', 'bold');
    
    fileLabel = uilabel(leftGrid, 'Text', 'No file loaded', ...
        'FontSize', 11, 'FontColor', [0.5 0.5 0.5]);
    
    infoLabel = uilabel(leftGrid, 'Text', '', 'FontSize', 10);
    
    uilabel(leftGrid, 'Text', '──── Data Processing ────', ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        'FontColor', [0.3 0.3 0.3]);
    
    cleanDataBtn = uibutton(leftGrid, 'Text', '⚕ Process Clinical Data', ...
        'ButtonPushedFcn', @(btn,event) cleanDataFunc(), ...
        'Enable', 'off', ...
        'BackgroundColor', [0.2 0.7 0.3], 'FontColor', 'white', ...
        'FontSize', 11, 'FontWeight', 'bold');
    
    cleanInfoLabel = uilabel(leftGrid, 'Text', '', 'FontSize', 10);
    
    uilabel(leftGrid, 'Text', '──── Patient Analysis ────', ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        'FontColor', [0.3 0.3 0.3]);
    
    uilabel(leftGrid, 'Text', 'Patient ID:', ...
        'FontWeight', 'bold');
    patientDropdown = uidropdown(leftGrid, 'Items', {'No data'}, ...
        'ValueChangedFcn', @(dd,event) updateDisplay(), ...
        'Enable', 'off', ...
        'Editable', 'on', ...
        'Placeholder', 'Type to search ID');
    
    uilabel(leftGrid, 'Text', 'Multiple IDs (comma-separated):', ...
        'FontSize', 9, 'FontColor', [0.4 0.4 0.4]);
    multiPatientField = uieditfield(leftGrid, 'text', ...
        'Placeholder', 'e.g., 2025050,2025051,2026007', ...
        'Enable', 'off', ...
        'ValueChangedFcn', @(src, event) updateDisplayMultiPatient());
    
    uilabel(leftGrid, 'Text', 'Sample Location:', ...
        'FontWeight', 'bold');
    catheterPlaceDropdown = uidropdown(leftGrid, 'Items', {'All Places'}, ...
        'Enable', 'off');
    
    uilabel(leftGrid, 'Text', 'Blood Gas Parameter:', ...
        'FontWeight', 'bold');
    paramDropdown = uidropdown(leftGrid, 'Items', paramNames, ...
        'Enable', 'off');
    
    plotBtn = uibutton(leftGrid, 'Text', '📈 Trend Analysis', ...
        'ButtonPushedFcn', @(btn,event) plotData(), ...
        'Enable', 'off', ...
        'BackgroundColor', [0.25 0.55 0.85], 'FontColor', 'white', ...
        'FontSize', 10, 'FontWeight', 'bold');
    
    multiPlotBtn = uibutton(leftGrid, 'Text', '📊 Multi-Parameter View', ...
        'ButtonPushedFcn', @(btn,event) plotMultiParameters(), ...
        'Enable', 'off', ...
        'BackgroundColor', [0.5 0.3 0.8], 'FontColor', 'white', ...
        'FontSize', 10, 'FontWeight', 'bold');
    
    placesPlotBtn = uibutton(leftGrid, 'Text', '🔬 Location Comparison', ...
        'ButtonPushedFcn', @(btn,event) plotParameterAcrossPlaces(), ...
        'Enable', 'off', ...
        'BackgroundColor', [0.9 0.5 0.2], 'FontColor', 'white', ...
        'FontSize', 10, 'FontWeight', 'bold');
    
    uilabel(leftGrid, 'Text', '──── Data Export ────', ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    exportCleanBtn = uibutton(leftGrid, 'Text', '💾 Export Cleaned Data', ...
        'ButtonPushedFcn', @(btn,event) exportCleanedData(), ...
        'Enable', 'off');
    
    exportPatientBtn = uibutton(leftGrid, 'Text', '📄 Export Patient Data', ...
        'ButtonPushedFcn', @(btn,event) exportPatientData(), ...
        'Enable', 'off');
    
    uilabel(leftGrid, 'Text', '─── Statistics ───', ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    statsBtn = uibutton(leftGrid, 'Text', '📊 Statistical Analysis', ...
        'ButtonPushedFcn', @(btn,event) calculateStatistics(), ...
        'Enable', 'off', ...
        'BackgroundColor', [0.2 0.6 0.8], 'FontColor', 'white', ...
        'FontSize', 10, 'FontWeight', 'bold');
    
    statusLabel = uilabel(leftGrid, 'Text', 'System Ready', ...
        'FontSize', 10, 'FontColor', [0 0.5 0], 'FontWeight', 'bold');
    
    %% RIGHT PANEL - Visualization
    rightPanel = uipanel(mainGrid, 'Title', 'Clinical Visualization & Data', ...
        'FontWeight', 'bold', 'FontSize', 11);
    rightGrid = uigridlayout(rightPanel, [2 1]);
    rightGrid.RowHeight = {'1x', '1x'};
    
    plotPanel = uipanel(rightGrid, 'Title', 'Graphical Analysis');
    ax = uiaxes(plotPanel, 'Position', [40 40 1100 320]);
    xlabel(ax, '');
    ylabel(ax, '');
    title(ax, 'Ready for data visualization');
    grid(ax, 'on');
    
    tablePanel = uipanel(rightGrid, 'Title', 'Patient Data Records');
    dataTable = uitable(tablePanel, 'Position', [10 10 1125 320]);
    
    %% Callback Functions
    
    function selectFile()
        [file, path] = uigetfile({'*.xlsx;*.csv', 'Excel and CSV Files (*.xlsx, *.csv)'}, 'Select PatLog File');
        if isequal(file, 0)
            return;
        end
        
        currentFile = fullfile(path, file);
        fileLabel.Text = sprintf('File: %s', file);
        statusLabel.Text = 'Loading file...';
        statusLabel.FontColor = [1 0.5 0];
        drawnow;
        
        try
            [~, ~, ext] = fileparts(currentFile);
            
            if strcmpi(ext, '.csv')
                fid = fopen(currentFile, 'r', 'n', 'UTF-8');
                
                if fid == -1
                    error('Could not open CSV file');
                end
                
                try
                    firstLine = fgetl(fid);
                    fclose(fid);
                    
                    numSemicolons = length(strfind(firstLine, ';'));
                    numCommas = length(strfind(firstLine, ','));
                    
                    if numSemicolons > numCommas
                        delimiter = ';';
                        numCols = numSemicolons + 1;
                        statusLabel.Text = 'Detected semicolon-delimited CSV';
                    else
                        delimiter = ',';
                        numCols = numCommas + 1;
                        statusLabel.Text = 'Detected comma-delimited CSV';
                    end
                    
                    drawnow;
                    
                    % Create format string for textscan
                    formatSpec = repmat('%q', 1, numCols);
                    
                    fid = fopen(currentFile, 'r', 'n', 'UTF-8');
                    dataCell = textscan(fid, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);
                    fclose(fid);
                    
                    numRows = length(dataCell{1});
                    rawData = cell(numRows, numCols);
                    
                    for col = 1:numCols
                        for row = 1:numRows
                            if col <= length(dataCell)
                                if row <= length(dataCell{col})
                                    val = dataCell{col}{row};
                                    if isempty(val)
                                        rawData{row, col} = '';
                                    else
                                        rawData{row, col} = val;
                                    end
                                else
                                    rawData{row, col} = '';
                                end
                            else
                                rawData{row, col} = '';
                            end
                        end
                    end
                    
                    rawData = [cell(1, size(rawData, 2)); rawData];
                    
                catch ME
                    if fid ~= -1
                        fclose(fid);
                    end
                    rethrow(ME);
                end
                
            else
                rawDataTemp = readcell(currentFile, 'Sheet', 1);
                
                rawData = rawDataTemp;
                for i = 1:numel(rawData)
                    if ismissing(rawData{i})
                        rawData{i} = '';
                    end
                end
            end
            
            infoLabel.Text = sprintf('Loaded: %d rows × %d cols', ...
                size(rawData, 1), size(rawData, 2));
            
            cleanDataBtn.Enable = 'on';
            
            statusLabel.Text = 'File loaded successfully';
            statusLabel.FontColor = [0 0.5 0];
            
            if size(rawData, 2) < 100
                statusLabel.Text = 'File appears pre-cleaned. Ready to analyze.';
                cleanData = rawData;
                enableAnalysis();
            end
            
        catch ME
            statusLabel.Text = sprintf('Error: %s', ME.message);
            statusLabel.FontColor = [1 0 0];
        end
    end
    
    function cleanDataFunc()
        if isempty(rawData)
            statusLabel.Text = 'No data to clean';
            statusLabel.FontColor = [1 0 0];
            return;
        end
        
        statusLabel.Text = 'Cleaning data...';
        statusLabel.FontColor = [1 0.5 0];
        drawnow;
        
        try
            columnsToKeep = {'Time', 'Sample #', 'Patient Id', 'Last Name', 'First Name', ...
                           'Sample type', 'Patient sex', '"pH"', 'pO2', 'pCO2', '"T"', ...
                           'tHb', 'sO2', 'Lac', 'K+', 'Na+', 'Cl-', 'FIO2'};
            
            headers = rawData(2, :);
            
            keepColIndices = [];
            foundColumns = {};
            
            for i = 1:length(columnsToKeep)
                colName = columnsToKeep{i};
                found = false;
                
                if strcmp(colName, '"T"')
                    if size(rawData, 2) >= 26
                        keepColIndices(end+1) = 26;
                        foundColumns{end+1} = char(string(headers{26}));
                        fprintf('Added Temperature at column 26: "%s"\n', char(string(headers{26})));
                        found = true;
                    else
                        fprintf('Warning: Column 26 (Temperature) does not exist\n');
                    end
                    continue;  % Skip to next column
                end
                
                for j = 1:length(headers)
                    if ~isempty(headers{j})
                        headerStr = strtrim(char(string(headers{j})));
                        
                        % Special handling for pH - must match exactly "pH" with quotes
                        if strcmp(colName, '"pH"')
                            if strcmp(headerStr, '"pH"') || strcmp(headerStr, 'pH')
                                if ~ismember(j, keepColIndices)
                                    keepColIndices(end+1) = j;
                                    foundColumns{end+1} = headerStr;
                                    fprintf('Found pH at column %d: "%s"\n', j, headerStr);
                                    found = true;
                                    break;
                                end
                            end
                        elseif strcmp(colName, 'pO2')
                            if contains(headerStr, 'pO2', 'IgnoreCase', true)
                                if ~contains(headerStr, 'pO2(v)', 'IgnoreCase', true)
                                    if ~ismember(j, keepColIndices)
                                        keepColIndices(end+1) = j;
                                        foundColumns{end+1} = headerStr;
                                        fprintf('Found pO2 at column %d: "%s"\n', j, headerStr);
                                        found = true;
                                        break;
                                    end
                                end
                            end
                        elseif strcmp(colName, 'sO2')
                            if contains(headerStr, 'sO2', 'IgnoreCase', true)
                                if ~contains(headerStr, 'sO2(v)', 'IgnoreCase', true)
                                    if ~ismember(j, keepColIndices)
                                        keepColIndices(end+1) = j;
                                        foundColumns{end+1} = headerStr;
                                        fprintf('Found sO2 at column %d: "%s"\n', j, headerStr);
                                        found = true;
                                        break;
                                    end
                                end
                            end
                        else
                            if contains(headerStr, colName, 'IgnoreCase', true)
                                if ~ismember(j, keepColIndices)
                                    keepColIndices(end+1) = j;
                                    foundColumns{end+1} = headerStr;
                                    found = true;
                                    break;
                                end
                            end
                        end
                    end
                end
                
                if ~found
                    fprintf('Warning: Column "%s" not found\n', colName);
                end
            end
            
            [keepColIndices, uniqueIdx] = unique(keepColIndices, 'stable');
            foundColumns = foundColumns(uniqueIdx);
            
            % Create cleaned data with only selected columns
            cleanData = rawData(:, keepColIndices);
            
            originalCols = size(rawData, 2);
            keptCols = length(keepColIndices);
            removedCols = originalCols - keptCols;
            
            cleanInfoLabel.Text = sprintf('Kept %d essential columns\nRemoved %d unnecessary columns\n\nColumns kept:\n%s', ...
                keptCols, removedCols, strjoin(foundColumns, ', '));
            
            statusLabel.Text = 'Data cleaned successfully';
            statusLabel.FontColor = [0 0.5 0];
            
            enableAnalysis();
            
        catch ME
            statusLabel.Text = sprintf('Cleaning error: %s', ME.message);
            statusLabel.FontColor = [1 0 0];
        end
    end
    
    function enableAnalysis()
        if isempty(cleanData)
            return;
        end
        
        headers = cleanData(2, :);
        
        fprintf('\n=== ENABLE ANALYSIS - Column Mapping ===\n');
        fprintf('Total columns after cleaning: %d\n', length(headers));
        fprintf('Headers:\n');
        for i = 1:length(headers)
            fprintf('  Col %d: "%s"\n', i, string(headers{i}));
        end
        
        paramIndices = zeros(1, length(paramNames));
        fprintf('\nSearching for parameters:\n');
        for i = 1:length(paramNames)
            % Special handling for Temperature
            if strcmp(paramNames{i}, 'T (°C)')
                for j = 1:length(headers)
                    if ~isempty(headers{j})
                        try
                            headerStr = char(string(headers{j}));
                            cleanHeader = strtrim(headerStr);
                            
                            if length(cleanHeader) >= 2 && cleanHeader(1) == '"' && cleanHeader(end) == '"'
                                cleanHeader = cleanHeader(2:end-1);
                            end
                            
                            if length(cleanHeader) >= 5 && ...
                               cleanHeader(1) == 'T' && ...
                               cleanHeader(2) == ' ' && ...
                               cleanHeader(3) == '(' && ...
                               contains(cleanHeader, 'C)')
                                paramIndices(i) = j;
                                fprintf('  ✓ %s found at column %d (%s)\n', paramNames{i}, j, headerStr);
                                break;
                            end
                        catch
                            continue;
                        end
                    end
                end
            else
                for j = 1:length(headers)
                    headerIsEmpty = isempty(headers{j});
                    if ~headerIsEmpty
                        try
                            headerStr = string(headers{j});
                            paramStr = string(paramNames{i});
                            if contains(headerStr, paramStr, 'IgnoreCase', true)
                                paramIndices(i) = j;
                                fprintf('  ✓ %s found at column %d (%s)\n', paramNames{i}, j, char(headerStr));
                                break;
                            end
                        catch
                            continue;
                        end
                    end
                end
            end
            if paramIndices(i) == 0
                fprintf('  ✗ %s NOT FOUND\n', paramNames{i});
            end
        end
        
        fprintf('\nParameter indices: %s\n', mat2str(paramIndices));
        
        dataRows = cleanData(3:end, :);
        for i = 1:length(paramNames)
            if paramIndices(i) > 0
                fprintf('  %s (col %d): ', paramNames{i}, paramIndices(i));
                for row = 1:min(3, size(dataRows, 1))
                    val = dataRows{row, paramIndices(i)};
                    fprintf('[%s] ', string(val));
                end
                fprintf('\n');
            end
        end
        fprintf('======================================\n\n');
        
        patientIDCol = 0;
        for j = 1:length(headers)
            if ~isempty(headers{j})
                try
                    headerStr = strtrim(string(headers{j}));
                    if strcmpi(headerStr, 'Patient Id') || strcmpi(headerStr, 'PatientId') || ...
                       strcmpi(headerStr, 'Patient ID') || strcmpi(headerStr, 'PatientID')
                        patientIDCol = j;
                        break;
                    end
                catch
                end
            end
        end
        
        if patientIDCol > 0
            dataRows = cleanData(3:end, :);
            patientIDs = [];
            
            for i = 1:size(dataRows, 1)
                val = dataRows{i, patientIDCol};
                if isnumeric(val)
                    if ~isnan(val)
                        if val ~= 0
                            patientIDs(end+1) = val;
                        end
                    end
                elseif ischar(val) || isstring(val)
                    if ~isempty(val)
                        numVal = str2double(val);
                        if ~isnan(numVal)
                            patientIDs(end+1) = numVal;
                        end
                    end
                end
            end
            
            if ~isempty(patientIDs)
                uniquePatients = unique(patientIDs);
                patientList = arrayfun(@num2str, uniquePatients, 'UniformOutput', false);
                patientDropdown.Items = ['All Patients', patientList];
                patientDropdown.Value = 'All Patients';
            else
                patientDropdown.Items = {'All Patients'};
            end
        else
            patientDropdown.Items = {'All Patients'};
        end
        
        catheterPlaceCol = 0;
        for i = 1:length(headers)
            if ~isempty(headers{i})
                try
                    headerStr = string(headers{i});
                    if contains(headerStr, 'First Name', 'IgnoreCase', true)
                        catheterPlaceCol = i;
                        break;
                    end
                catch
                end
            end
        end
        
        if catheterPlaceCol > 0
            dataRows = cleanData(3:end, :);
            catheterPlaces = {};
            for i = 1:size(dataRows, 1)
                val = dataRows{i, catheterPlaceCol};
                if ischar(val) || isstring(val)
                    valStr = strtrim(char(val));
                    if ~isempty(valStr)
                        catheterPlaces{end+1} = valStr;
                    end
                end
            end
            uniquePlaces = unique(catheterPlaces);
            catheterPlaceDropdown.Items = ['All Places', uniquePlaces];
            catheterPlaceDropdown.Value = 'All Places';
        end
        
        foundParamNames = {};
        for i = 1:length(paramNames)
            if paramIndices(i) > 0
                foundParamNames{end+1} = paramNames{i};
            end
        end
        if ~isempty(foundParamNames)
            paramDropdown.Items = foundParamNames;
            paramDropdown.Value = foundParamNames{1};
        end
        
        patientDropdown.Enable = 'on';
        multiPatientField.Enable = 'on';
        catheterPlaceDropdown.Enable = 'on';
        paramDropdown.Enable = 'on';
        plotBtn.Enable = 'on';
        multiPlotBtn.Enable = 'on';
        placesPlotBtn.Enable = 'on';
        exportCleanBtn.Enable = 'on';
        exportPatientBtn.Enable = 'on';
        statsBtn.Enable = 'on';
        
        updateDisplay();
    end
    
    function updateDisplay()
        if isempty(cleanData)
            return;
        end
        
        selectedPatient = patientDropdown.Value;
        dataRows = cleanData(3:end, :);
        
        % Filter by patient
        keepAllRows = strcmp(selectedPatient, 'All Patients');
        
        if ~keepAllRows
            if patientIDCol > 0
                patientNum = str2double(selectedPatient);
                keepRows = false(size(dataRows, 1), 1);
                for i = 1:size(dataRows, 1)
                    val = dataRows{i, patientIDCol};
                    numVal = toNumber(val);
                    if ~isnan(numVal)
                        if numVal == patientNum
                            keepRows(i) = true;
                        end
                    end
                end
                dataRows = dataRows(keepRows, :);
                dataRows = sortRowsByTime(dataRows, headers);
                displayData = [cleanData(1:2, :); dataRows];
            else
                dataRows = sortRowsByTime(dataRows, headers);
                displayData = [cleanData(1:2, :); dataRows];
            end
        else
            dataRows = sortRowsByTime(dataRows, headers);
            displayData = [cleanData(1:2, :); dataRows];
        end
        
        updateCatheterPlaceDropdown(dataRows);
        
        % Show all rows in the table
        tableData = displayData;
        
        for i = 1:numel(tableData)
            try
                if ismissing(tableData{i})
                    tableData{i} = '';
                end
            catch
            end
        end
        
        dataTable.Data = tableData;
        dataTable.ColumnName = 'numbered';
        
        statusLabel.Text = sprintf('Displaying %d rows (filtered)', size(displayData, 1) - 2);
        statusLabel.FontColor = [0 0.5 0];
    end
    
    function updateDisplayMultiPatient()
        if isempty(cleanData)
            return;
        end
        
        multiIDText = strtrim(multiPatientField.Value);
        
        if isempty(multiIDText)
            % If empty, revert to single patient dropdown
            updateDisplay();
            return;
        end
        
        % Parse comma-separated IDs
        idStrings = strsplit(multiIDText, ',');
        patientNums = [];
        for i = 1:length(idStrings)
            idStr = strtrim(idStrings{i});
            num = str2double(idStr);
            if ~isnan(num)
                patientNums(end+1) = num;
            end
        end
        
        if isempty(patientNums)
            statusLabel.Text = 'Invalid patient IDs';
            statusLabel.FontColor = [1 0 0];
            return;
        end
        
        dataRows = cleanData(3:end, :);
        
        % Filter by multiple patients
        keepRows = false(size(dataRows, 1), 1);
        for i = 1:size(dataRows, 1)
            val = dataRows{i, patientIDCol};
            numVal = toNumber(val);
            if ~isnan(numVal)
                for p = 1:length(patientNums)
                    if numVal == patientNums(p)
                        keepRows(i) = true;
                        break;
                    end
                end
            end
        end
        
        dataRows = dataRows(keepRows, :);
        dataRows = sortRowsByTime(dataRows, headers);
        displayData = [cleanData(1:2, :); dataRows];
        
        updateCatheterPlaceDropdown(dataRows);
        
        % Show all rows in the table
        tableData = displayData;
        
        for i = 1:numel(tableData)
            try
                if ismissing(tableData{i})
                    tableData{i} = '';
                end
            catch
            end
        end
        
        dataTable.Data = tableData;
        dataTable.ColumnName = 'numbered';
        
        statusLabel.Text = sprintf('Displaying %d rows from %d patients', ...
            size(displayData, 1) - 2, length(patientNums));
        statusLabel.FontColor = [0 0.5 0];
    end
    
    function updateCatheterPlaceDropdown(dataRows)
        catheterPlaceCol = 0;
        for i = 1:length(headers)
            if ~isempty(headers{i})
                try
                    headerStr = string(headers{i});
                    if contains(headerStr, 'First Name', 'IgnoreCase', true)
                        catheterPlaceCol = i;
                        break;
                    end
                catch
                end
            end
        end
        
        if catheterPlaceCol > 0
            catheterPlaces = {};
            for i = 1:size(dataRows, 1)
                val = dataRows{i, catheterPlaceCol};
                if ischar(val) || isstring(val)
                    valStr = strtrim(char(val));
                    if ~isempty(valStr)
                        catheterPlaces{end+1} = valStr;
                    end
                end
            end
            uniquePlaces = unique(catheterPlaces);
            
            currentSelection = catheterPlaceDropdown.Value;
            
            catheterPlaceDropdown.Items = ['All Places', uniquePlaces];
            
            if ismember(currentSelection, catheterPlaceDropdown.Items)
                catheterPlaceDropdown.Value = currentSelection;
            else
                catheterPlaceDropdown.Value = 'All Places';
            end
        end
    end
    
    function plotData()
        if isempty(cleanData)
            statusLabel.Text = 'No data loaded';
            statusLabel.FontColor = [1 0 0];
            return;
        end
        
        % Get selected patients (prioritize multi-field)
        [patientNums, patientLabel] = getSelectedPatients();
        
        if isempty(patientNums) && ~strcmp(patientLabel, 'All Patients')
            statusLabel.Text = 'No valid patient IDs selected';
            statusLabel.FontColor = [1 0 0];
            return;
        end
        
        selectedParam = paramDropdown.Value;
        selectedCatheterPlace = catheterPlaceDropdown.Value;
        
        % Check if multiple patients selected
        if length(patientNums) > 1
            % Create separate window with subplots for each patient
            plotMultiPatientTrend(patientNums, selectedParam, selectedCatheterPlace);
        else
            % Single patient - plot in main axes
            plotSinglePatientTrend(patientLabel, selectedParam, selectedCatheterPlace, ax);
        end
    end
    
    function plotMultiPatientTrend(patientNums, selectedParam, selectedCatheterPlace)
        % Create new figure with subplots for multiple patients
        numPatients = length(patientNums);
        rows = ceil(sqrt(numPatients));
        cols = ceil(numPatients / rows);
        
        paramWithUnit = getParamWithUnit(selectedParam);
        
        trendFig = figure('Name', sprintf('%s Trend - Multiple Patients', paramWithUnit), ...
            'Position', [50 50 1400 800]);
        
        for p = 1:numPatients
            subplot(rows, cols, p);
            plotSinglePatientTrendInAxes(num2str(patientNums(p)), selectedParam, selectedCatheterPlace, gca);
        end
        
        sgtitle(sprintf('%s Across Patients - Location: %s', paramWithUnit, selectedCatheterPlace), ...
            'FontSize', 16, 'FontWeight', 'bold');
        
        statusLabel.Text = sprintf('Plotted %s for %d patients', selectedParam, numPatients);
        statusLabel.FontColor = [0 0.5 0];
    end
    
    function plotSinglePatientTrend(patientLabel, selectedParam, selectedCatheterPlace, axes_handle)
        plotSinglePatientTrendInAxes(patientLabel, selectedParam, selectedCatheterPlace, axes_handle);
        statusLabel.Text = sprintf('Plotted %s for patient %s', selectedParam, patientLabel);
        statusLabel.FontColor = [0 0.5 0];
    end
    
    function plotSinglePatientTrendInAxes(patientLabel, selectedParam, selectedCatheterPlace, axes_handle)
        % Core plotting logic for a single patient
        paramIdx = 0;
        for i = 1:length(paramNames)
            if strcmp(paramNames{i}, selectedParam)
                paramIdx = i;
                break;
            end
        end
        
        if paramIdx == 0 || paramIndices(paramIdx) == 0
            title(axes_handle, sprintf('Patient %s - Parameter not available', patientLabel));
            return;
        end
        
        colIdx = paramIndices(paramIdx);
        dataRows = cleanData(3:end, :);
        
        % Filter by patient
        if ~strcmp(patientLabel, 'All Patients')
            if patientIDCol > 0
                patientNum = str2double(patientLabel);
                keepRows = false(size(dataRows, 1), 1);
                for i = 1:size(dataRows, 1)
                    val = dataRows{i, patientIDCol};
                    numVal = toNumber(val);
                    if ~isnan(numVal) && numVal == patientNum
                        keepRows(i) = true;
                    end
                end
                dataRows = dataRows(keepRows, :);
            end
        end
        
        % Filter by catheter place
        if ~strcmp(selectedCatheterPlace, 'All Places')
            catheterPlaceCol = 0;
            for i = 1:length(headers)
                if ~isempty(headers{i})
                    try
                        headerStr = string(headers{i});
                        if contains(headerStr, 'First Name', 'IgnoreCase', true)
                            catheterPlaceCol = i;
                            break;
                        end
                    catch
                    end
                end
            end
            
            if catheterPlaceCol > 0
                keepRows = false(size(dataRows, 1), 1);
                for i = 1:size(dataRows, 1)
                    val = dataRows{i, catheterPlaceCol};
                    if ischar(val) || isstring(val)
                        if strcmp(strtrim(char(val)), selectedCatheterPlace)
                            keepRows(i) = true;
                        end
                    end
                end
                dataRows = dataRows(keepRows, :);
            end
        end
        
        % Sort chronologically so plot reads left-to-right in time
        dataRows = sortRowsByTime(dataRows, headers);
        
        % Get time column
        timeCol = 0;
        for i = 1:length(headers)
            if ~isempty(headers{i})
                try
                    headerStr = string(headers{i});
                    if contains(headerStr, 'Last Name', 'IgnoreCase', true)
                        timeCol = i;
                        break;
                    end
                catch
                end
            end
        end
        
        % Extract values
        values = [];
        indices = [];
        
        for i = 1:size(dataRows, 1)
            val = dataRows{i, colIdx};
            numVal = toNumber(val);
            
            if ~isnan(numVal)
                values(end+1) = numVal;
                indices(end+1) = i;
            end
        end
        
        % Plot
        if isempty(values)
            title(axes_handle, sprintf('Patient %s - No data', patientLabel), ...
                'FontWeight', 'bold');
            text(axes_handle, 0.5, 0.5, 'No data available', ...
                'HorizontalAlignment', 'center', 'FontSize', 14);
            set(axes_handle, 'XTick', [], 'YTick', []);
        else
            plot(axes_handle, indices, values, '-o', 'LineWidth', 2.5, 'MarkerSize', 8, ...
                'MarkerFaceColor', [0.2 0.5 0.8], 'Color', [0.2 0.5 0.8]);
            
            % Use shared helper for x-axis labels
            [tickPts, tickLbls] = buildTimeAxisLabels(dataRows, headers);
            if ~isempty(tickPts)
                xticks(axes_handle, tickPts);
                xticklabels(axes_handle, tickLbls);
                xtickangle(axes_handle, 45);
            end
            xlabel(axes_handle, 'Time', 'FontWeight', 'bold');
            
            paramWithUnit = getParamWithUnit(selectedParam);
            ylabel(axes_handle, paramWithUnit, 'FontWeight', 'bold');
            title(axes_handle, sprintf('Patient %s - %s', patientLabel, paramWithUnit), ...
                'FontWeight', 'bold');
            grid(axes_handle, 'on');
        end
    end
    
    function [tickPoints, tickLabels] = buildTimeAxisLabels(dataRows, hdrs)
        % Build x-axis tick labels for location/trend plots:
        % - Named ticks from Last Name (first occurrence of each unique label)
        % - Between named ticks, show HH:MM from Time column
        
        % Find Last Name column
        lastNameCol = 0;
        for i = 1:length(hdrs)
            if ~isempty(hdrs{i})
                try
                    hStr = string(hdrs{i});
                    if contains(hStr, 'Last Name', 'IgnoreCase', true)
                        lastNameCol = i;
                        break;
                    end
                catch
                end
            end
        end
        
        % Find actual Time column
        actTimeCol = 0;
        for i = 1:length(hdrs)
            if ~isempty(hdrs{i})
                try
                    hStr = string(hdrs{i});
                    if strcmpi(strtrim(hStr), 'Time')
                        actTimeCol = i;
                        break;
                    end
                catch
                end
            end
        end
        
        tickPoints = [];
        tickLabels = {};
        seenLabels = containers.Map('KeyType', 'char', 'ValueType', 'double');
        
        for i = 1:size(dataRows, 1)
            lastNameVal = '';
            timeHHMM = '';
            
            % Check Last Name column
            if lastNameCol > 0
                tv = dataRows{i, lastNameCol};
                if ischar(tv) && ~isempty(strtrim(tv))
                    lastNameVal = strtrim(tv);
                elseif isstring(tv) && strlength(strtrim(tv)) > 0
                    lastNameVal = char(strtrim(tv));
                end
            end
            
            % Get HH:MM from actual Time column
            if actTimeCol > 0
                tVal = dataRows{i, actTimeCol};
                if isnumeric(tVal) && ~isnan(tVal)
                    try
                        timeHHMM = datestr(tVal, 'HH:MM');
                    catch
                        timeHHMM = num2str(tVal);
                    end
                elseif ischar(tVal) && ~isempty(strtrim(tVal))
                    rawTime = strtrim(tVal);
                    try
                        dt = datetime(rawTime);
                        timeHHMM = datestr(dt, 'HH:MM');
                    catch
                        tokens = regexp(rawTime, '(\d{1,2}:\d{2})', 'tokens');
                        if ~isempty(tokens)
                            timeHHMM = tokens{1}{1};
                        else
                            timeHHMM = rawTime;
                        end
                    end
                elseif isstring(tVal) && strlength(strtrim(tVal)) > 0
                    rawTime = char(strtrim(tVal));
                    try
                        dt = datetime(rawTime);
                        timeHHMM = datestr(dt, 'HH:MM');
                    catch
                        tokens = regexp(rawTime, '(\d{1,2}:\d{2})', 'tokens');
                        if ~isempty(tokens)
                            timeHHMM = tokens{1}{1};
                        else
                            timeHHMM = rawTime;
                        end
                    end
                elseif isdatetime(tVal) && ~isnat(tVal)
                    timeHHMM = datestr(tVal, 'HH:MM');
                end
            end
            
            % Decide which label to use
            if ~isempty(lastNameVal)
                if ~isKey(seenLabels, lastNameVal)
                    seenLabels(lastNameVal) = i;
                    tickPoints(end+1) = i;
                    tickLabels{end+1} = lastNameVal;
                end
            elseif ~isempty(timeHHMM)
                if ~isKey(seenLabels, timeHHMM)
                    seenLabels(timeHHMM) = i;
                    tickPoints(end+1) = i;
                    tickLabels{end+1} = timeHHMM;
                end
            end
        end
    end
    
    function sortedRows = sortRowsByTime(dataRows, hdrs)
        % Sort rows ascending by the actual Time column so plots read
        % left-to-right in chronological order. Rows whose Time can't be
        % parsed are placed at the end in their original order (stable sort).
        if isempty(dataRows)
            sortedRows = dataRows;
            return;
        end
        
        % Find actual Time column (same logic as buildTimeAxisLabels)
        actTimeCol = 0;
        for i = 1:length(hdrs)
            if ~isempty(hdrs{i})
                try
                    hStr = string(hdrs{i});
                    if strcmpi(strtrim(hStr), 'Time')
                        actTimeCol = i;
                        break;
                    end
                catch
                end
            end
        end
        
        if actTimeCol == 0
            sortedRows = dataRows;
            return;
        end
        
        n = size(dataRows, 1);
        sortKey = NaN(n, 1);
        
        for i = 1:n
            tVal = dataRows{i, actTimeCol};
            dt = NaT;
            try
                if isnumeric(tVal) && ~isnan(tVal)
                    dt = datetime(tVal, 'ConvertFrom', 'datenum');
                elseif isdatetime(tVal) && ~isnat(tVal)
                    dt = tVal;
                elseif (ischar(tVal) || isstring(tVal)) && strlength(string(tVal)) > 0
                    raw = strtrim(char(tVal));
                    parsed = false;
                    % Try common Czech locale format first: DD.MM.YYYY HH:MM
                    fmts = {'dd.MM.yyyy HH:mm', 'dd.MM.yyyy HH:mm:ss', ...
                            'dd/MM/yyyy HH:mm', 'yyyy-MM-dd HH:mm:ss', ...
                            'yyyy-MM-dd HH:mm', 'MM/dd/yyyy HH:mm'};
                    for f = 1:length(fmts)
                        try
                            dt = datetime(raw, 'InputFormat', fmts{f});
                            parsed = true;
                            break;
                        catch
                        end
                    end
                    if ~parsed
                        try
                            dt = datetime(raw);
                        catch
                            dt = NaT;
                        end
                    end
                end
            catch
                dt = NaT;
            end
            
            if ~isnat(dt)
                sortKey(i) = datenum(dt);
            end
        end
        
        % Stable ascending sort: NaNs (unparseable) go to the end
        [~, order] = sortrows([isnan(sortKey), sortKey], [1 2]);
        sortedRows = dataRows(order, :);
    end
    
    function exportCleanedData()
        if isempty(cleanData)
            statusLabel.Text = 'No cleaned data to export';
            statusLabel.FontColor = [1 0 0];
            return;
        end
        
        [file, path] = uiputfile('*.xlsx', 'Save Cleaned Data As', 'PatLog_cleaned.xlsx');
        if isequal(file, 0)
            return;
        end
        
        try
            fullPath = fullfile(path, file);
            writecell(cleanData, fullPath);
            statusLabel.Text = sprintf('Exported to: %s', file);
            statusLabel.FontColor = [0 0.5 0];
        catch ME
            statusLabel.Text = sprintf('Export error: %s', ME.message);
            statusLabel.FontColor = [1 0 0];
        end
    end
    
    function exportPatientData()
        if isempty(cleanData)
            statusLabel.Text = 'No data to export';
            statusLabel.FontColor = [1 0 0];
            return;
        end
        
        selectedPatient = patientDropdown.Value;
        
        if strcmp(selectedPatient, 'All Patients')
            statusLabel.Text = 'Select a specific patient to export';
            statusLabel.FontColor = [1 0 0];
            return;
        end
        
        patientNum = str2double(selectedPatient);
        dataRows = cleanData(3:end, :);
        
        fprintf('\n=== EXPORT PATIENT DATA DEBUG ===\n');
        fprintf('Selected patient: %s (numeric: %.0f)\n', selectedPatient, patientNum);
        fprintf('Patient ID column: %d\n', patientIDCol);
        
        % Filter data
        keepRows = false(size(dataRows, 1), 1);
        matchCount = 0;
        for i = 1:size(dataRows, 1)
            val = dataRows{i, patientIDCol};
            numVal = toNumber(val);
            if i <= 3
                fprintf('Row %d: Patient ID = "%s" (class: %s) -> numeric: %.0f\n', ...
                    i, char(string(val)), class(val), numVal);
            end
            if ~isnan(numVal)
                if numVal == patientNum
                    keepRows(i) = true;
                    matchCount = matchCount + 1;
                    if matchCount <= 3
                    end
                end
            end
        end
        
        fprintf('Total matching rows: %d\n', matchCount);
        fprintf('================================\n');
        
        exportData = [cleanData(1:2, :); dataRows(keepRows, :)];
        
        [file, path] = uiputfile('*.xlsx', 'Save Patient Data As', ...
            sprintf('PatLog_patient_%s.xlsx', selectedPatient));
        if isequal(file, 0)
            return;
        end
        
        try
            fullPath = fullfile(path, file);
            writecell(exportData, fullPath);
            statusLabel.Text = sprintf('Exported patient %s to: %s', selectedPatient, file);
            statusLabel.FontColor = [0 0.5 0];
        catch ME
            statusLabel.Text = sprintf('Export error: %s', ME.message);
            statusLabel.FontColor = [1 0 0];
        end
    end
    
    function plotMultiParameters()
        if isempty(cleanData)
            return;
        end
        
        % Get selected patients (prioritize multi-field)
        [patientNums, patientLabel] = getSelectedPatients();
        
        % Determine if multiple patients are selected
        multiPatientMode = length(patientNums) > 1;
        
        % Find catheter place column to get available locations
        cathPlaceCol = 0;
        for i = 1:length(headers)
            if ~isempty(headers{i})
                try
                    headerStr = string(headers{i});
                    if contains(headerStr, 'First Name', 'IgnoreCase', true)
                        cathPlaceCol = i;
                        break;
                    end
                catch
                end
            end
        end
        
        % Get unique locations per patient and determine shared locations
        availableLocations = {};
        locationPatientMap = containers.Map('KeyType', 'char', 'ValueType', 'any');  % location -> list of patient nums
        sharedLocations = {};  % locations present in ALL selected patients
        
        if cathPlaceCol > 0
            dataRows = cleanData(3:end, :);
            
            if ~isempty(patientNums) && patientIDCol > 0
                % Build per-patient location sets
                for pIdx = 1:length(patientNums)
                    for i = 1:size(dataRows, 1)
                        pidVal = dataRows{i, patientIDCol};
                        pidNum = toNumber(pidVal);
                        if ~isnan(pidNum) && pidNum == patientNums(pIdx)
                            locVal = dataRows{i, cathPlaceCol};
                            locStr = '';
                            if ischar(locVal) && ~isempty(strtrim(locVal))
                                locStr = strtrim(locVal);
                            elseif isstring(locVal) && strlength(strtrim(locVal)) > 0
                                locStr = char(strtrim(locVal));
                            end
                            if ~isempty(locStr)
                                if isKey(locationPatientMap, locStr)
                                    existing = locationPatientMap(locStr);
                                    if ~ismember(patientNums(pIdx), existing)
                                        locationPatientMap(locStr) = [existing, patientNums(pIdx)];
                                    end
                                else
                                    locationPatientMap(locStr) = patientNums(pIdx);
                                end
                            end
                        end
                    end
                end
                availableLocations = sort(keys(locationPatientMap));
                
                % Find shared locations (present in ALL selected patients)
                if length(patientNums) > 1
                    for lIdx = 1:length(availableLocations)
                        loc = availableLocations{lIdx};
                        pList = locationPatientMap(loc);
                        if length(pList) == length(patientNums)
                            sharedLocations{end+1} = loc;
                        end
                    end
                end
            else
                % No patient filtering - get all locations
                for i = 1:size(dataRows, 1)
                    locVal = dataRows{i, cathPlaceCol};
                    if ischar(locVal) && ~isempty(strtrim(locVal))
                        availableLocations{end+1} = strtrim(locVal);
                    elseif isstring(locVal) && strlength(strtrim(locVal)) > 0
                        availableLocations{end+1} = char(strtrim(locVal));
                    end
                end
                availableLocations = unique(availableLocations);
            end
        end
        
        % Calculate grid dimensions
        numLocations = length(availableLocations);
        numParams = length(paramNames);
        numPatients = 0;
        if multiPatientMode
            numPatients = length(patientNums);
        end
        
        % Layout: 4 columns
        % Row 1: Patient title (col 1-2) | Location title (col 3-4)
        % Rows 2-N: Patient checkboxes (col 1-2) | Location checkboxes (col 3-4)
        % Next row: separator across all columns
        % Next row: Parameter title across all columns
        % Next rows: Parameters in 4-column grid
        % Last rows: Plot + Cancel buttons
        
        locRowsNeeded = ceil((numLocations + 1) / 2);  % +1 for "All Places"
        patRowsNeeded = max(1, ceil(numPatients / 2));
        leftRightRows = max(patRowsNeeded, locRowsNeeded);
        paramGridRows = ceil(numParams / 4);
        hasLegend = multiPatientMode && ~isempty(sharedLocations) && length(sharedLocations) < numLocations;
        
        hasOverlapRows = multiPatientMode * 2;  % two rows for overlap checkboxes
        totalGridRows = 1 + leftRightRows + hasOverlapRows + hasLegend + 1 + 1 + paramGridRows + 2;
        numCols = 4;
        
        paramSelectionFig = uifigure('Name', 'Select Patients, Locations and Parameters', ...
            'Position', [150 80 800 min(750, max(400, 35 * totalGridRows))]);
        
        g = uigridlayout(paramSelectionFig, [totalGridRows, numCols]);
        g.ColumnWidth = {'1x', '1x', '1x', '1x'};
        g.RowHeight = repmat({'fit'}, 1, totalGridRows);
        g.Padding = [10 10 10 10];
        g.RowSpacing = 4;
        g.ColumnSpacing = 8;
        g.Scrollable = 'on';
        
        % Define colors for location coding
        sharedColor = [0.0 0.5 0.0];    % Green for shared across all patients
        individualColor = [0.6 0.3 0.0]; % Brown/orange for individual patient only
        
        currentRow = 1;
        
        % === Row 1: Section titles ===
        if multiPatientMode
            patTitle = uilabel(g, 'Text', 'Select patients:', ...
                'FontWeight', 'bold', 'FontSize', 13, 'FontColor', [0.2 0.4 0.7]);
            patTitle.Layout.Row = currentRow;
            patTitle.Layout.Column = [1 2];
        else
            patTitle = uilabel(g, 'Text', '', 'FontSize', 1);
            patTitle.Layout.Row = currentRow;
            patTitle.Layout.Column = [1 2];
        end
        
        locTitle = uilabel(g, 'Text', 'Select sample location(s):', ...
            'FontWeight', 'bold', 'FontSize', 13, 'FontColor', [0.2 0.4 0.7]);
        locTitle.Layout.Row = currentRow;
        locTitle.Layout.Column = [3 4];
        
        currentRow = currentRow + 1;
        
        % === Rows 2-N: Patient checkboxes (left) | Location checkboxes (right) ===
        patientCheckboxes = {};
        if multiPatientMode
            patientCheckboxes = cell(1, numPatients);
            for i = 1:numPatients
                r = currentRow + floor((i-1) / 2);
                c = mod(i-1, 2) + 1;
                patientCheckboxes{i} = uicheckbox(g, ...
                    'Text', sprintf('ID: %d', patientNums(i)), ...
                    'Value', true);
                patientCheckboxes{i}.Layout.Row = r;
                patientCheckboxes{i}.Layout.Column = c;
            end
        end
        
        % Location checkboxes on right side (columns 3-4)
        % First: "All Places" checkbox
        allPlacesCheckbox = uicheckbox(g, 'Text', 'All Places', 'Value', true, ...
            'FontWeight', 'bold', ...
            'ValueChangedFcn', @(cb, event) toggleAllPlaces(cb.Value));
        allPlacesCheckbox.Layout.Row = currentRow;
        allPlacesCheckbox.Layout.Column = 3;
        
        locationCheckboxes = cell(1, numLocations);
        for i = 1:numLocations
            r = currentRow + floor(i / 2);
            c = mod(i, 2) + 3;  % columns 3-4
            
            % Determine if this location is shared or individual
            isShared = false;
            if multiPatientMode
                for sIdx = 1:length(sharedLocations)
                    if strcmp(availableLocations{i}, sharedLocations{sIdx})
                        isShared = true;
                        break;
                    end
                end
            end
            
            if multiPatientMode && ~isempty(sharedLocations) && length(sharedLocations) < numLocations
                if isShared
                    labelColor = sharedColor;
                else
                    labelColor = individualColor;
                end
            else
                labelColor = [0 0 0];  % default black
            end
            
            locationCheckboxes{i} = uicheckbox(g, ...
                'Text', availableLocations{i}, ...
                'Value', false, ...
                'FontColor', labelColor, ...
                'ValueChangedFcn', @(cb, event) onLocationChecked(cb.Value));
            locationCheckboxes{i}.Layout.Row = r;
            locationCheckboxes{i}.Layout.Column = c;
        end
        
        currentRow = currentRow + leftRightRows;
        
        % === Overlap Trends checkboxes (only in multi-patient mode) ===
        overlapCheckbox = [];
        simpleOverlapCheckbox = [];
        if multiPatientMode
            simpleOverlapCheckbox = uicheckbox(g, ...
                'Text', 'Simple overlap (same graph, independent x-axes)', ...
                'Value', false, ...
                'FontWeight', 'bold', 'FontColor', [0.6 0.0 0.0]);
            simpleOverlapCheckbox.Layout.Row = currentRow;
            simpleOverlapCheckbox.Layout.Column = [1 2];
            currentRow = currentRow + 1;
            
            overlapCheckbox = uicheckbox(g, ...
                'Text', 'Aligned overlap (compare on shared timeline)', ...
                'Value', false, ...
                'FontWeight', 'bold', 'FontColor', [0.0 0.0 0.6]);
            overlapCheckbox.Layout.Row = currentRow;
            overlapCheckbox.Layout.Column = [1 2];
            currentRow = currentRow + 1;
            
            % Make checkboxes mutually exclusive
            simpleOverlapCheckbox.ValueChangedFcn = @(cb, event) onSimpleOverlapChanged(cb.Value);
            overlapCheckbox.ValueChangedFcn = @(cb, event) onAlignedOverlapChanged(cb.Value);
        end
        
        % === Legend row (if applicable) ===
        if hasLegend
            sharedLegend = uilabel(g, 'Text', '● Shared by all patients', ...
                'FontSize', 10, 'FontColor', sharedColor, 'FontWeight', 'bold');
            sharedLegend.Layout.Row = currentRow;
            sharedLegend.Layout.Column = 3;
            
            indivLegend = uilabel(g, 'Text', '● Individual patient only', ...
                'FontSize', 10, 'FontColor', individualColor, 'FontWeight', 'bold');
            indivLegend.Layout.Row = currentRow;
            indivLegend.Layout.Column = 4;
            
            currentRow = currentRow + 1;
        end
        
        % === Separator row ===
        sep = uilabel(g, 'Text', ' ', 'FontSize', 2);
        sep.Layout.Row = currentRow;
        sep.Layout.Column = [1 4];
        currentRow = currentRow + 1;
        
        % === Parameter title ===
        paramTitle = uilabel(g, 'Text', 'Select parameters to plot:', ...
            'FontWeight', 'bold', 'FontSize', 13, 'FontColor', [0.2 0.4 0.7]);
        paramTitle.Layout.Row = currentRow;
        paramTitle.Layout.Column = [1 4];
        currentRow = currentRow + 1;
        
        % === Parameters in 4-column grid ===
        paramCheckboxes = cell(1, numParams);
        plotParams = {};  % No parameters pre-selected
        
        for i = 1:numParams
            r = currentRow + floor((i-1) / 4);
            c = mod(i-1, 4) + 1;
            isChecked = false;
            for j = 1:length(plotParams)
                if strcmp(paramNames{i}, plotParams{j})
                    isChecked = true;
                    break;
                end
            end
            paramCheckboxes{i} = uicheckbox(g, 'Text', paramNames{i}, 'Value', isChecked);
            paramCheckboxes{i}.Layout.Row = r;
            paramCheckboxes{i}.Layout.Column = c;
        end
        
        currentRow = currentRow + paramGridRows;
        
        % === Buttons ===
        plotMultiBtn = uibutton(g, 'Text', 'Plot Selected', ...
            'ButtonPushedFcn', @(btn,event) executeMultiPlot(), ...
            'BackgroundColor', [0.2 0.8 0.4], 'FontColor', 'white', ...
            'FontWeight', 'bold');
        plotMultiBtn.Layout.Row = currentRow;
        plotMultiBtn.Layout.Column = [1 2];
        
        cancelBtn = uibutton(g, 'Text', 'Cancel', ...
            'ButtonPushedFcn', @(btn,event) close(paramSelectionFig));
        cancelBtn.Layout.Row = currentRow;
        cancelBtn.Layout.Column = [3 4];
        
        function onLocationChecked(isChecked)
            if isChecked
                allPlacesCheckbox.Value = false;
            end
        end
        
        function onSimpleOverlapChanged(isChecked)
            if isChecked && ~isempty(overlapCheckbox)
                overlapCheckbox.Value = false;
            end
        end
        
        function onAlignedOverlapChanged(isChecked)
            if isChecked && ~isempty(simpleOverlapCheckbox)
                simpleOverlapCheckbox.Value = false;
            end
        end
        
        function toggleAllPlaces(isAllSelected)
            % When "All Places" is checked, uncheck individual locations
            % When unchecked, leave individual selections as they are
            if isAllSelected
                for lIdx = 1:length(locationCheckboxes)
                    locationCheckboxes{lIdx}.Value = false;
                end
            end
        end
        
        function executeMultiPlot()
            % Get selected parameters
            selectedParams = {};
            selectedIndices = [];
            for k = 1:length(paramCheckboxes)
                if paramCheckboxes{k}.Value
                    if paramIndices(k) > 0
                        selectedParams{end+1} = paramNames{k};
                        selectedIndices(end+1) = paramIndices(k);
                    end
                end
            end
            
            if isempty(selectedParams)
                uialert(paramSelectionFig, 'Please select at least one parameter', 'No Selection');
                return;
            end
            
            % Get selected locations
            selectedLocations = {};
            useAllPlaces = allPlacesCheckbox.Value;
            if ~useAllPlaces
                for lIdx = 1:length(locationCheckboxes)
                    if locationCheckboxes{lIdx}.Value
                        selectedLocations{end+1} = availableLocations{lIdx};
                    end
                end
                if isempty(selectedLocations)
                    % Nothing selected individually and All Places unchecked
                    useAllPlaces = true;
                end
            end
            
            if useAllPlaces
                selectedCatheterPlace = 'All Places';
            elseif length(selectedLocations) == 1
                selectedCatheterPlace = selectedLocations{1};
            else
                selectedCatheterPlace = selectedLocations;  % cell array of multiple
            end
            
            % Get selected patients (only in multi-patient mode)
            selectedPatientNums = [];
            if multiPatientMode
                for k = 1:length(patientCheckboxes)
                    if patientCheckboxes{k}.Value
                        selectedPatientNums(end+1) = patientNums(k);
                    end
                end
                
                if isempty(selectedPatientNums)
                    uialert(paramSelectionFig, 'Please select at least one patient', 'No Selection');
                    return;
                end
            else
                % Single patient mode
                if strcmp(patientLabel, 'All Patients')
                    selectedPatientNums = [];  % All patients
                else
                    selectedPatientNums = str2double(patientLabel);
                end
            end
            
            % Check overlap modes BEFORE closing figure (UI objects destroyed on close)
            overlapMode = false;
            simpleOverlapMode = false;
            if ~isempty(overlapCheckbox) && overlapCheckbox.Value
                overlapMode = true;
            end
            if ~isempty(simpleOverlapCheckbox) && simpleOverlapCheckbox.Value
                simpleOverlapMode = true;
            end
            
            close(paramSelectionFig);
            
            % Execute the multi-parameter plot
            plotMultiParamsForPatients(selectedPatientNums, selectedParams, selectedIndices, selectedCatheterPlace, patientLabel, overlapMode, simpleOverlapMode);
        end
    end
    
    function plotMultiParamsForPatients(patientNums, selectedParams, selectedIndices, selectedCatheterPlace, patientLabel, overlapMode, simpleOverlapMode)
        if nargin < 6
            overlapMode = false;
        end
        if nargin < 7
            simpleOverlapMode = false;
        end
        
        % Check if multiple patients
        if length(patientNums) > 1
            
            if simpleOverlapMode
                % ============================================================
                % SIMPLE OVERLAP: Each patient plotted 1..N on same axes
                % No timeline alignment - just visual comparison of trends
                % ============================================================
                numPatients = length(patientNums);
                numParams = length(selectedParams);
                
                % Define distinct colors for each patient
                patientColors = [
                    0.0  0.45 0.74;   % Blue
                    0.85 0.33 0.10;   % Red-orange
                    0.47 0.67 0.19;   % Green
                    0.49 0.18 0.56;   % Purple
                    0.93 0.69 0.13;   % Gold
                    0.30 0.75 0.93;   % Cyan
                    0.64 0.08 0.18;   % Dark red
                    0.00 0.62 0.45;   % Teal
                ];
                if numPatients > size(patientColors, 1)
                    patientColors = [patientColors; lines(numPatients - size(patientColors, 1))];
                end
                
                % Marker styles for extra differentiation
                markerStyles = {'o', 's', 'd', '^', 'v', '>', '<', 'p'};
                
                % Collect filtered data per patient
                patientDataSets = cell(1, numPatients);
                patientTickPts = cell(1, numPatients);
                patientTickLbls = cell(1, numPatients);
                
                catheterPlaceCol = 0;
                for i = 1:length(headers)
                    if ~isempty(headers{i})
                        try
                            headerStr = string(headers{i});
                            if contains(headerStr, 'First Name', 'IgnoreCase', true)
                                catheterPlaceCol = i;
                                break;
                            end
                        catch
                        end
                    end
                end
                
                for p = 1:numPatients
                    dataRows = cleanData(3:end, :);
                    
                    % Filter by patient
                    keepRows = false(size(dataRows, 1), 1);
                    for i = 1:size(dataRows, 1)
                        val = dataRows{i, patientIDCol};
                        numVal = toNumber(val);
                        if ~isnan(numVal) && numVal == patientNums(p)
                            keepRows(i) = true;
                        end
                    end
                    dataRows = dataRows(keepRows, :);
                    
                    % Filter by location
                    skipFilter = false;
                    if ischar(selectedCatheterPlace) && strcmp(selectedCatheterPlace, 'All Places')
                        skipFilter = true;
                    end
                    
                    if ~skipFilter && catheterPlaceCol > 0
                        keepRows = false(size(dataRows, 1), 1);
                        for i = 1:size(dataRows, 1)
                            val = dataRows{i, catheterPlaceCol};
                            if ischar(val) || isstring(val)
                                valStr = strtrim(char(val));
                                if iscell(selectedCatheterPlace)
                                    for locIdx = 1:length(selectedCatheterPlace)
                                        if strcmp(valStr, selectedCatheterPlace{locIdx})
                                            keepRows(i) = true;
                                            break;
                                        end
                                    end
                                else
                                    if strcmp(valStr, selectedCatheterPlace)
                                        keepRows(i) = true;
                                    end
                                end
                            end
                        end
                        dataRows = dataRows(keepRows, :);
                    end
                    
                    patientDataSets{p} = dataRows;
                    % Sort chronologically so plot reads left-to-right in time
                    patientDataSets{p} = sortRowsByTime(patientDataSets{p}, headers);
                    [tp, tl] = buildTimeAxisLabels(patientDataSets{p}, headers);
                    patientTickPts{p} = tp;
                    patientTickLbls{p} = tl;
                end
                
                % Create figure
                rows = ceil(sqrt(numParams));
                cols = ceil(numParams / rows);
                
                patientInfo = sprintf('Patients: %s (Simple Overlap)', ...
                    strjoin(arrayfun(@num2str, patientNums, 'UniformOutput', false), ', '));
                
                multiParamFig = figure('Name', sprintf('Multi-Parameter Analysis - %s', patientInfo), ...
                    'Position', [30 30 min(1800, max(1000, 450*cols)), min(1050, max(650, 350*rows))]);
                
                for k = 1:numParams
                    ax = subplot(rows, cols, k);
                    hold(ax, 'on');
                    colIdx = selectedIndices(k);
                    legendEntries = {};
                    legendHandles = [];
                    
                    for p = 1:numPatients
                        dataRows = patientDataSets{p};
                        
                        values = [];
                        indices = [];
                        for i = 1:size(dataRows, 1)
                            val = dataRows{i, colIdx};
                            numVal = toNumber(val);
                            if ~isnan(numVal)
                                values(end+1) = numVal;
                                indices(end+1) = i;
                            end
                        end
                        
                        if ~isempty(values)
                            mrkr = markerStyles{mod(p-1, length(markerStyles)) + 1};
                            h = plot(ax, indices, values, ['-' mrkr], ...
                                'LineWidth', 2.5, 'MarkerSize', 7, ...
                                'Color', patientColors(p, :), ...
                                'MarkerFaceColor', patientColors(p, :));
                            legendHandles(end+1) = h;
                            legendEntries{end+1} = sprintf('ID: %d', patientNums(p));
                        end
                    end
                    
                    hold(ax, 'off');
                    
                    % Build combined x-axis labels colored by patient
                    % Collect all tick positions and labels across patients
                    % Each patient's labels are drawn as colored text
                    set(ax, 'XTickLabel', []);
                    
                    % Collect all unique tick positions
                    allTickPts = [];
                    for p = 1:numPatients
                        allTickPts = union(allTickPts, patientTickPts{p});
                    end
                    if ~isempty(allTickPts)
                        xticks(ax, allTickPts);
                    end
                    
                    % Draw colored labels
                    yLim = ylim(ax);
                    yRange = yLim(2) - yLim(1);
                    
                    % Stack labels for different patients at slightly different y offsets
                    for p = 1:numPatients
                        yOff = yLim(1) - (0.02 + (p-1) * 0.05) * yRange;
                        tp = patientTickPts{p};
                        tl = patientTickLbls{p};
                        for ti = 1:length(tp)
                            text(ax, tp(ti), yOff, tl{ti}, ...
                                'HorizontalAlignment', 'right', ...
                                'VerticalAlignment', 'top', ...
                                'Rotation', 45, ...
                                'FontSize', 7, 'FontWeight', 'bold', ...
                                'Color', patientColors(p, :), ...
                                'Clipping', 'off');
                        end
                    end
                    
                    paramWithUnit = getParamWithUnit(selectedParams{k});
                    ylabel(ax, paramWithUnit, 'FontSize', 10, 'FontWeight', 'bold');
                    title(ax, paramWithUnit, 'FontSize', 11, 'FontWeight', 'bold');
                    grid(ax, 'on');
                    
                    if ~isempty(legendHandles)
                        legend(ax, legendHandles, legendEntries, ...
                            'Location', 'best', 'FontSize', 8);
                    end
                end
                
                if iscell(selectedCatheterPlace)
                    locLabel = strjoin(selectedCatheterPlace, ', ');
                else
                    locLabel = selectedCatheterPlace;
                end
                sgtitle(sprintf('Simple Overlap - Location: %s', locLabel), ...
                    'FontSize', 13, 'FontWeight', 'bold');
                
                % Adjust subplot positions for room
                allAxes = findobj(multiParamFig, 'Type', 'axes');
                for aIdx = 1:length(allAxes)
                    pos = allAxes(aIdx).Position;
                    newBottom = pos(2) + 0.08;
                    newHeight = pos(4) - 0.10;
                    allAxes(aIdx).Position = [pos(1), newBottom, pos(3), newHeight];
                end
                
            elseif overlapMode
                % ============================================================
                % OVERLAPPED MODE: All patients on same axes per parameter
                % ============================================================
                numPatients = length(patientNums);
                numParams = length(selectedParams);
                
                % Define distinct colors for each patient
                patientColors = [
                    0.0  0.45 0.74;   % Blue
                    0.85 0.33 0.10;   % Red-orange
                    0.47 0.67 0.19;   % Green
                    0.49 0.18 0.56;   % Purple
                    0.93 0.69 0.13;   % Gold
                    0.30 0.75 0.93;   % Cyan
                    0.64 0.08 0.18;   % Dark red
                    0.00 0.62 0.45;   % Teal
                ];
                if numPatients > size(patientColors, 1)
                    patientColors = [patientColors; lines(numPatients - size(patientColors, 1))];
                end
                
                % Step 1: Collect data per patient (filtered by location)
                patientDataSets = cell(1, numPatients);
                patientLastNames = cell(1, numPatients);  % Last Name labels per patient
                
                catheterPlaceCol = 0;
                for i = 1:length(headers)
                    if ~isempty(headers{i})
                        try
                            headerStr = string(headers{i});
                            if contains(headerStr, 'First Name', 'IgnoreCase', true)
                                catheterPlaceCol = i;
                                break;
                            end
                        catch
                        end
                    end
                end
                
                lastNameCol = 0;
                for i = 1:length(headers)
                    if ~isempty(headers{i})
                        try
                            headerStr = string(headers{i});
                            if contains(headerStr, 'Last Name', 'IgnoreCase', true)
                                lastNameCol = i;
                                break;
                            end
                        catch
                        end
                    end
                end
                
                actTimeCol = 0;
                for i = 1:length(headers)
                    if ~isempty(headers{i})
                        try
                            headerStr = string(headers{i});
                            if strcmpi(strtrim(headerStr), 'Time')
                                actTimeCol = i;
                                break;
                            end
                        catch
                        end
                    end
                end
                
                for p = 1:numPatients
                    dataRows = cleanData(3:end, :);
                    
                    % Filter by patient
                    keepRows = false(size(dataRows, 1), 1);
                    for i = 1:size(dataRows, 1)
                        val = dataRows{i, patientIDCol};
                        numVal = toNumber(val);
                        if ~isnan(numVal) && numVal == patientNums(p)
                            keepRows(i) = true;
                        end
                    end
                    dataRows = dataRows(keepRows, :);
                    
                    % Filter by location
                    skipFilter = false;
                    if ischar(selectedCatheterPlace) && strcmp(selectedCatheterPlace, 'All Places')
                        skipFilter = true;
                    end
                    
                    if ~skipFilter && catheterPlaceCol > 0
                        keepRows = false(size(dataRows, 1), 1);
                        for i = 1:size(dataRows, 1)
                            val = dataRows{i, catheterPlaceCol};
                            if ischar(val) || isstring(val)
                                valStr = strtrim(char(val));
                                if iscell(selectedCatheterPlace)
                                    for locIdx = 1:length(selectedCatheterPlace)
                                        if strcmp(valStr, selectedCatheterPlace{locIdx})
                                            keepRows(i) = true;
                                            break;
                                        end
                                    end
                                else
                                    if strcmp(valStr, selectedCatheterPlace)
                                        keepRows(i) = true;
                                    end
                                end
                            end
                        end
                        dataRows = dataRows(keepRows, :);
                    end
                    
                    patientDataSets{p} = dataRows;
                    
                    % Extract Last Name labels for this patient
                    lastNames = {};
                    if lastNameCol > 0
                        for i = 1:size(dataRows, 1)
                            tv = dataRows{i, lastNameCol};
                            lnStr = '';
                            if ischar(tv) && ~isempty(strtrim(tv))
                                lnStr = strtrim(tv);
                            elseif isstring(tv) && strlength(strtrim(tv)) > 0
                                lnStr = char(strtrim(tv));
                            end
                            lastNames{end+1} = lnStr;
                        end
                    end
                    patientLastNames{p} = lastNames;
                end
                
                % Step 2: Build unified timeline from all patients' Last Name labels
                % Collect all unique Last Name labels in order of first appearance
                unifiedLabels = {};  % ordered unique labels
                labelPatientSource = {};  % which patient(s) each label came from
                
                % Merge labels from all patients while preserving relative order
                % Use a simple approach: iterate through each patient's labels and insert new ones
                for p = 1:numPatients
                    lNames = patientLastNames{p};
                    for i = 1:length(lNames)
                        if ~isempty(lNames{i})
                            alreadyExists = false;
                            for u = 1:length(unifiedLabels)
                                if strcmp(unifiedLabels{u}, lNames{i})
                                    alreadyExists = true;
                                    % Add this patient as a source
                                    if ~ismember(p, labelPatientSource{u})
                                        labelPatientSource{u}(end+1) = p;
                                    end
                                    break;
                                end
                            end
                            if ~alreadyExists
                                unifiedLabels{end+1} = lNames{i};
                                labelPatientSource{end+1} = p;
                            end
                        end
                    end
                end
                
                numUnifiedLabels = length(unifiedLabels);
                
                % Also add HH:MM time points between named labels for each patient
                % For overlap, we map each patient's data row to the nearest unified label index
                % Rows with named labels map exactly; rows without get intermediate positions
                
                % Build per-patient x-position mapping
                patientXPositions = cell(1, numPatients);  % x-coords for each patient's data rows
                
                for p = 1:numPatients
                    dataRows = patientDataSets{p};
                    nRows = size(dataRows, 1);
                    xPos = nan(1, nRows);
                    
                    lNames = patientLastNames{p};
                    
                    % Find positions of named labels in this patient's data
                    namedRowIdx = [];  % row indices that have a Last Name
                    namedUnifiedIdx = [];  % corresponding unified label index
                    
                    for i = 1:nRows
                        if i <= length(lNames) && ~isempty(lNames{i})
                            % Find this label in unified labels
                            for u = 1:numUnifiedLabels
                                if strcmp(unifiedLabels{u}, lNames{i})
                                    namedRowIdx(end+1) = i;
                                    namedUnifiedIdx(end+1) = u;
                                    xPos(i) = u;
                                    break;
                                end
                            end
                        end
                    end
                    
                    % Interpolate positions for rows between named labels
                    if length(namedRowIdx) >= 2
                        for seg = 1:(length(namedRowIdx)-1)
                            startRow = namedRowIdx(seg);
                            endRow = namedRowIdx(seg+1);
                            startX = namedUnifiedIdx(seg);
                            endX = namedUnifiedIdx(seg+1);
                            
                            gapRows = endRow - startRow;
                            if gapRows > 1
                                for ri = (startRow+1):(endRow-1)
                                    fraction = (ri - startRow) / gapRows;
                                    xPos(ri) = startX + fraction * (endX - startX);
                                end
                            end
                        end
                    elseif length(namedRowIdx) == 1
                        % Only one named label - offset others around it
                        anchorRow = namedRowIdx(1);
                        anchorX = namedUnifiedIdx(1);
                        for i = 1:nRows
                            if isnan(xPos(i))
                                xPos(i) = anchorX + (i - anchorRow) * 0.3;
                            end
                        end
                    end
                    
                    % Handle rows before first named label or after last
                    if ~isempty(namedRowIdx)
                        firstNamedRow = namedRowIdx(1);
                        firstNamedX = namedUnifiedIdx(1);
                        for i = 1:(firstNamedRow-1)
                            if isnan(xPos(i))
                                xPos(i) = firstNamedX - (firstNamedRow - i) * 0.5;
                            end
                        end
                        
                        lastNamedRow = namedRowIdx(end);
                        lastNamedX = namedUnifiedIdx(end);
                        for i = (lastNamedRow+1):nRows
                            if isnan(xPos(i))
                                xPos(i) = lastNamedX + (i - lastNamedRow) * 0.5;
                            end
                        end
                    else
                        % No named labels at all - just use sequential numbering
                        for i = 1:nRows
                            xPos(i) = i;
                        end
                    end
                    
                    patientXPositions{p} = xPos;
                end
                
                % Step 3: Create the overlapped figure
                rows = ceil(sqrt(numParams));
                cols = ceil(numParams / rows);
                
                patientInfo = sprintf('Patients: %s (Overlapped)', ...
                    strjoin(arrayfun(@num2str, patientNums, 'UniformOutput', false), ', '));
                
                multiParamFig = figure('Name', sprintf('Multi-Parameter Analysis - %s', patientInfo), ...
                    'Position', [30 30 min(1800, max(1000, 450*cols)), min(1050, max(650, 350*rows))]);
                
                for k = 1:numParams
                    ax = subplot(rows, cols, k);
                    hold(ax, 'on');
                    colIdx = selectedIndices(k);
                    legendEntries = {};
                    legendHandles = [];
                    
                    for p = 1:numPatients
                        dataRows = patientDataSets{p};
                        xPos = patientXPositions{p};
                        
                        values = [];
                        xCoords = [];
                        for i = 1:size(dataRows, 1)
                            val = dataRows{i, colIdx};
                            numVal = toNumber(val);
                            if ~isnan(numVal) && ~isnan(xPos(i))
                                values(end+1) = numVal;
                                xCoords(end+1) = xPos(i);
                            end
                        end
                        
                        if ~isempty(values)
                            h = plot(ax, xCoords, values, '-o', ...
                                'LineWidth', 2.5, 'MarkerSize', 7, ...
                                'Color', patientColors(p, :), ...
                                'MarkerFaceColor', patientColors(p, :));
                            legendHandles(end+1) = h;
                            legendEntries{end+1} = sprintf('ID: %d', patientNums(p));
                        end
                    end
                    
                    hold(ax, 'off');
                    
                    % Set unified x-axis tick labels
                    if ~isempty(unifiedLabels)
                        xticks(ax, 1:numUnifiedLabels);
                        xlim(ax, [0.5, numUnifiedLabels + 0.5]);
                        
                        % Color each x-tick label by patient source
                        % Clear default labels and draw colored text manually
                        set(ax, 'XTickLabel', []);
                        
                        yLim = ylim(ax);
                        yRange = yLim(2) - yLim(1);
                        % Place labels just below the axis bottom edge
                        yTextPos = yLim(1) - 0.02 * yRange;
                        
                        for u = 1:numUnifiedLabels
                            % Determine color: shared = black, single patient = patient color
                            if length(labelPatientSource{u}) == 1
                                txtColor = patientColors(labelPatientSource{u}(1), :);
                            else
                                txtColor = [0 0 0];  % black for shared labels
                            end
                            
                            text(ax, u, yTextPos, unifiedLabels{u}, ...
                                'HorizontalAlignment', 'right', ...
                                'VerticalAlignment', 'top', ...
                                'Rotation', 45, ...
                                'FontSize', 7, 'FontWeight', 'bold', ...
                                'Color', txtColor, ...
                                'Clipping', 'off');
                        end
                    end
                    
                    paramWithUnit = getParamWithUnit(selectedParams{k});
                    ylabel(ax, paramWithUnit, 'FontSize', 10, 'FontWeight', 'bold');
                    title(ax, paramWithUnit, 'FontSize', 11, 'FontWeight', 'bold');
                    grid(ax, 'on');
                    
                    if ~isempty(legendHandles)
                        legend(ax, legendHandles, legendEntries, ...
                            'Location', 'best', 'FontSize', 8);
                    end
                end
                
                if iscell(selectedCatheterPlace)
                    locLabel = strjoin(selectedCatheterPlace, ', ');
                else
                    locLabel = selectedCatheterPlace;
                end
                sgtitle(sprintf('Overlapped Multi-Parameter Analysis - Location: %s', locLabel), ...
                    'FontSize', 13, 'FontWeight', 'bold');
                
                % Adjust subplot positions to avoid title overlap and leave room for labels
                allAxes = findobj(multiParamFig, 'Type', 'axes');
                for aIdx = 1:length(allAxes)
                    pos = allAxes(aIdx).Position;
                    % Raise bottom, shrink height: more space at bottom for labels + key
                    newBottom = pos(2) + 0.08;
                    newTop = pos(4) - 0.10;
                    allAxes(aIdx).Position = [pos(1), newBottom, pos(3), newTop];
                end
                
                % Add color legend at bottom of figure using individual colored annotations
                totalItems = numPatients + 1;  % patients + "shared" label
                itemWidth = 0.11;
                prefixWidth = 0.07;
                sharedWidth = 0.19;
                totalWidth = prefixWidth + totalItems * itemWidth + sharedWidth;
                startX = max(0.02, (1.0 - totalWidth) / 2);
                
                annotation(multiParamFig, 'textbox', [startX 0.015 prefixWidth 0.03], ...
                    'String', 'X-axis color key:', ...
                    'HorizontalAlignment', 'left', ...
                    'FontSize', 9, 'FontWeight', 'bold', ...
                    'EdgeColor', 'none', 'FitBoxToText', 'off', ...
                    'Color', [0 0 0]);
                
                xOff = startX + prefixWidth;
                squareChar = char(9632);  % Unicode filled square
                for p = 1:numPatients
                    annotation(multiParamFig, 'textbox', [xOff 0.015 itemWidth 0.03], ...
                        'String', sprintf('%s ID:%d', squareChar, patientNums(p)), ...
                        'HorizontalAlignment', 'left', ...
                        'FontSize', 9, 'FontWeight', 'bold', ...
                        'EdgeColor', 'none', 'FitBoxToText', 'off', ...
                        'Color', patientColors(p, :));
                    xOff = xOff + itemWidth;
                end
                
                % Add shared label in black
                annotation(multiParamFig, 'textbox', [xOff 0.015 sharedWidth 0.03], ...
                    'String', sprintf('| %s Black = shared action', squareChar), ...
                    'HorizontalAlignment', 'left', ...
                    'FontSize', 9, 'FontWeight', 'bold', ...
                    'EdgeColor', 'none', 'FitBoxToText', 'off', ...
                    'Color', [0 0 0]);
            
            else
            numPatients = length(patientNums);
            numParams = length(selectedParams);
            
            % Calculate layout: patients as rows, parameters as columns
            rows = numPatients;
            cols = numParams;
            
            % Create title showing patient info
            patientInfo = sprintf('Patients: %s', strjoin(arrayfun(@num2str, patientNums, 'UniformOutput', false), ', '));
            
            multiParamFig = figure('Name', sprintf('Multi-Parameter Analysis - %s', patientInfo), ...
                'Position', [50 50 min(1600, 400*cols), min(1000, 300*rows)]);
            
            % Create subplots - one row per patient
            for p = 1:numPatients
                % Get data for this patient only
                dataRows = cleanData(3:end, :);
                
                % Filter by this patient
                keepRows = false(size(dataRows, 1), 1);
                for i = 1:size(dataRows, 1)
                    val = dataRows{i, patientIDCol};
                    numVal = toNumber(val);
                    if ~isnan(numVal) && numVal == patientNums(p)
                        keepRows(i) = true;
                    end
                end
                dataRows = dataRows(keepRows, :);
                
                % Filter by catheter place
                skipFilter = false;
                if ischar(selectedCatheterPlace) && strcmp(selectedCatheterPlace, 'All Places')
                    skipFilter = true;
                end
                
                if ~skipFilter
                    catheterPlaceCol = 0;
                    for i = 1:length(headers)
                        if ~isempty(headers{i})
                            try
                                headerStr = string(headers{i});
                                if contains(headerStr, 'First Name', 'IgnoreCase', true)
                                    catheterPlaceCol = i;
                                    break;
                                end
                            catch
                            end
                        end
                    end
                    
                    if catheterPlaceCol > 0
                        keepRows = false(size(dataRows, 1), 1);
                        for i = 1:size(dataRows, 1)
                            val = dataRows{i, catheterPlaceCol};
                            if ischar(val) || isstring(val)
                                valStr = strtrim(char(val));
                                if iscell(selectedCatheterPlace)
                                    % Multiple locations selected
                                    for locIdx = 1:length(selectedCatheterPlace)
                                        if strcmp(valStr, selectedCatheterPlace{locIdx})
                                            keepRows(i) = true;
                                            break;
                                        end
                                    end
                                else
                                    % Single location
                                    if strcmp(valStr, selectedCatheterPlace)
                                        keepRows(i) = true;
                                    end
                                end
                            end
                        end
                        dataRows = dataRows(keepRows, :);
                    end
                end
                
                % Sort chronologically so plot reads left-to-right in time
                dataRows = sortRowsByTime(dataRows, headers);
                
                % Get time column
                timeCol = 0;
                for i = 1:length(headers)
                    if ~isempty(headers{i})
                        try
                            headerStr = string(headers{i});
                            if contains(headerStr, 'Last Name', 'IgnoreCase', true)
                                timeCol = i;
                                break;
                            end
                        catch
                        end
                    end
                end
                
                % Build time axis labels using shared helper
                [tickPts, tickLbls] = buildTimeAxisLabels(dataRows, headers);
                
                % Plot parameters for this patient
                colors = lines(numParams);
                
                for k = 1:numParams
                    subplot(rows, cols, (p-1)*cols + k);
                    colIdx = selectedIndices(k);
                    
                    values = [];
                    indices = [];
                    for i = 1:size(dataRows, 1)
                        val = dataRows{i, colIdx};
                        numVal = toNumber(val);
                        
                        if ~isnan(numVal)
                            values(end+1) = numVal;
                            indices(end+1) = i;
                        end
                    end
                    
                    if ~isempty(values)
                        plot(indices, values, '-o', ...
                            'LineWidth', 2.5, 'MarkerSize', 8, ...
                            'Color', colors(k, :), ...
                            'MarkerFaceColor', colors(k, :));
                        
                        if ~isempty(tickPts)
                            xticks(tickPts);
                            xticklabels(tickLbls);
                            xtickangle(45);
                        end
                        
                        xlabel('Time', 'FontSize', 9);
                        paramWithUnit = getParamWithUnit(selectedParams{k});
                        ylabel(paramWithUnit, 'FontSize', 10, 'FontWeight', 'bold');
                        
                        % Title shows patient ID + parameter
                        title(sprintf('Patient %d - %s', patientNums(p), paramWithUnit), ...
                            'FontSize', 11, 'FontWeight', 'bold');
                        grid on;
                    else
                        title(sprintf('Patient %d - %s (No Data)', patientNums(p), selectedParams{k}), ...
                            'FontSize', 11);
                    end
                end
            end
            
            if iscell(selectedCatheterPlace)
                locLabel = strjoin(selectedCatheterPlace, ', ');
            else
                locLabel = selectedCatheterPlace;
            end
            sgtitle(sprintf('Multi-Parameter Analysis - Location: %s', locLabel), ...
                'FontSize', 14, 'FontWeight', 'bold');
            
            end  % end overlap if/else
            
        else
            % Single patient - original behavior (all params in one figure)
            dataRows = cleanData(3:end, :);
            
            % Filter by patient
            if ~isempty(patientNums)
                keepRows = false(size(dataRows, 1), 1);
                for i = 1:size(dataRows, 1)
                    val = dataRows{i, patientIDCol};
                    numVal = toNumber(val);
                    if ~isnan(numVal)
                        for p = 1:length(patientNums)
                            if numVal == patientNums(p)
                                keepRows(i) = true;
                                break;
                            end
                        end
                    end
                end
                dataRows = dataRows(keepRows, :);
            end
            
            % Filter by catheter place
            skipFilter = false;
            if ischar(selectedCatheterPlace) && strcmp(selectedCatheterPlace, 'All Places')
                skipFilter = true;
            end
            
            if ~skipFilter
                catheterPlaceCol = 0;
                for i = 1:length(headers)
                    if ~isempty(headers{i})
                        try
                            headerStr = string(headers{i});
                            if contains(headerStr, 'First Name', 'IgnoreCase', true)
                                catheterPlaceCol = i;
                                break;
                            end
                        catch
                        end
                    end
                end
                
                if catheterPlaceCol > 0
                    keepRows = false(size(dataRows, 1), 1);
                    for i = 1:size(dataRows, 1)
                        val = dataRows{i, catheterPlaceCol};
                        if ischar(val) || isstring(val)
                            valStr = strtrim(char(val));
                            if iscell(selectedCatheterPlace)
                                for locIdx = 1:length(selectedCatheterPlace)
                                    if strcmp(valStr, selectedCatheterPlace{locIdx})
                                        keepRows(i) = true;
                                        break;
                                    end
                                end
                            else
                                if strcmp(valStr, selectedCatheterPlace)
                                    keepRows(i) = true;
                                end
                            end
                        end
                    end
                    dataRows = dataRows(keepRows, :);
                end
            end
            
            % Sort chronologically so plot reads left-to-right in time
            dataRows = sortRowsByTime(dataRows, headers);
            
            % Get time column
            timeCol = 0;
            for i = 1:length(headers)
                if ~isempty(headers{i})
                    try
                        headerStr = string(headers{i});
                        if contains(headerStr, 'Last Name', 'IgnoreCase', true)
                            timeCol = i;
                            break;
                        end
                    catch
                    end
                end
            end
            
            % Build time axis labels using shared helper
            [tickPts, tickLbls] = buildTimeAxisLabels(dataRows, headers);
            
            % Create figure with subplots
            numParams = length(selectedParams);
            rows = ceil(sqrt(numParams));
            cols = ceil(numParams / rows);
            
            multiParamFig = figure('Name', sprintf('Multi-Parameter Analysis - Patient %s', patientLabel), ...
                'Position', [100 100 1400 900]);
            
            colors = lines(numParams);
            
            for k = 1:numParams
                subplot(rows, cols, k);
                colIdx = selectedIndices(k);
                
                values = [];
                indices = [];
                for i = 1:size(dataRows, 1)
                    val = dataRows{i, colIdx};
                    numVal = toNumber(val);
                    
                    if ~isnan(numVal)
                        values(end+1) = numVal;
                        indices(end+1) = i;
                    end
                end
                
                if ~isempty(values)
                    plot(indices, values, '-o', ...
                        'LineWidth', 2.5, 'MarkerSize', 8, ...
                        'Color', colors(k, :), ...
                        'MarkerFaceColor', colors(k, :));
                    
                    if ~isempty(tickPts)
                        xticks(tickPts);
                        xticklabels(tickLbls);
                        xtickangle(45);
                    end
                    
                    xlabel('Time', 'FontSize', 10);
                    paramWithUnit = getParamWithUnit(selectedParams{k});
                    ylabel(paramWithUnit, 'FontSize', 11, 'FontWeight', 'bold');
                    title(paramWithUnit, 'FontSize', 12, 'FontWeight', 'bold');
                    grid on;
                end
            end
            
            if iscell(selectedCatheterPlace)
                locLabel2 = strjoin(selectedCatheterPlace, ', ');
            else
                locLabel2 = selectedCatheterPlace;
            end
            sgtitle(sprintf('Patient %s - Location: %s', patientLabel, locLabel2), ...
                'FontSize', 14, 'FontWeight', 'bold');
        end
        
        statusLabel.Text = sprintf('Plotted %d parameters', length(selectedParams));
        statusLabel.FontColor = [0 0.5 0];
    end
    
    function plotParameterAcrossPlaces()
        if isempty(cleanData)
            return;
        end
        
        % Get selected patients (prioritize multi-field)
        [patientNums, patientLabel] = getSelectedPatients();
        
        if isempty(patientNums) && ~strcmp(patientLabel, 'All Patients')
            statusLabel.Text = 'No valid patient IDs selected';
            statusLabel.FontColor = [1 0 0];
            return;
        end
        
        selectedParam = paramDropdown.Value;
        
        % Check if multiple patients selected
        if length(patientNums) > 1
            % Show patient selection dialog for multi-patient mode
            showLocationComparisonMultiDialog(patientNums, selectedParam);
        else
            % Single patient mode - use original logic
            plotParameterAcrossPlacesSingle(patientLabel, selectedParam);
        end
    end
    
    function showLocationComparisonMultiDialog(patientNums, selectedParam)
        % Create dialog for patient selection
        dialogFig = uifigure('Name', 'Select Patients for Location Comparison', ...
            'Position', [300 200 500 400]);
        
        g = uigridlayout(dialogFig, [length(patientNums) + 4, 1]);
        g.RowHeight = repmat({'fit'}, 1, length(patientNums) + 4);
        
        uilabel(g, 'Text', sprintf('Select patients to compare %s across locations:', selectedParam), ...
            'FontWeight', 'bold', 'FontSize', 14, ...
            'FontColor', [0.2 0.4 0.7]);
        
        uilabel(g, 'Text', 'Each patient will be shown in a different color', ...
            'FontSize', 10, 'FontColor', [0.5 0.5 0.5]);
        
        % Patient checkboxes
        patientCheckboxes = cell(1, length(patientNums));
        for i = 1:length(patientNums)
            patientCheckboxes{i} = uicheckbox(g, ...
                'Text', sprintf('Patient ID: %d', patientNums(i)), ...
                'Value', true);
        end
        
        plotBtn = uibutton(g, 'Text', 'Compare Selected Patients', ...
            'ButtonPushedFcn', @(btn,event) executeLocationComparison(), ...
            'BackgroundColor', [0.2 0.8 0.4], 'FontColor', 'white', ...
            'FontWeight', 'bold');
        
        cancelBtn = uibutton(g, 'Text', 'Cancel', ...
            'ButtonPushedFcn', @(btn,event) close(dialogFig));
        
        function executeLocationComparison()
            % Get selected patients
            selectedPatientNums = [];
            for k = 1:length(patientCheckboxes)
                if patientCheckboxes{k}.Value
                    selectedPatientNums(end+1) = patientNums(k);
                end
            end
            
            if isempty(selectedPatientNums)
                uialert(dialogFig, 'Please select at least one patient', 'No Selection');
                return;
            end
            
            close(dialogFig);
            
            % Plot location comparison with multiple patients
            plotLocationComparisonMultiPatients(selectedPatientNums, selectedParam);
        end
    end
    
    function plotLocationComparisonMultiPatients(patientNums, selectedParam)
        % Multi-patient location comparison - separate subplot for each patient
        paramIdx = 0;
        for i = 1:length(paramNames)
            if strcmp(paramNames{i}, selectedParam)
                paramIdx = i;
                break;
            end
        end
        
        if paramIdx == 0 || paramIndices(paramIdx) == 0
            statusLabel.Text = 'Parameter not available';
            statusLabel.FontColor = [1 0 0];
            return;
        end
        
        colIdx = paramIndices(paramIdx);
        paramWithUnit = getParamWithUnit(selectedParam);
        
        % Find catheter place column (First Name)
        catheterPlaceCol = 0;
        for i = 1:length(headers)
            if ~isempty(headers{i})
                try
                    headerStr = string(headers{i});
                    if contains(headerStr, 'First Name', 'IgnoreCase', true)
                        catheterPlaceCol = i;
                        break;
                    end
                catch
                end
            end
        end
        
        if catheterPlaceCol == 0
            statusLabel.Text = 'First Name column not found';
            statusLabel.FontColor = [1 0 0];
            return;
        end
        
        % Find time column (Last Name)
        timeCol = 0;
        for i = 1:length(headers)
            if ~isempty(headers{i})
                try
                    headerStr = string(headers{i});
                    if contains(headerStr, 'Last Name', 'IgnoreCase', true)
                        timeCol = i;
                        break;
                    end
                catch
                end
            end
        end
        
        % Find actual Time column as fallback
        actualTimeCol = 0;
        for i = 1:length(headers)
            if ~isempty(headers{i})
                try
                    headerStr = string(headers{i});
                    if strcmpi(strtrim(headerStr), 'Time')
                        actualTimeCol = i;
                        break;
                    end
                catch
                end
            end
        end
        
        % Calculate layout for subplots
        numPatients = length(patientNums);
        rows = ceil(sqrt(numPatients));
        cols = ceil(numPatients / rows);
        
        % Create figure
        placeFig = figure('Name', sprintf('%s Across Locations - Multiple Patients', paramWithUnit), ...
            'Position', [50 50 min(1600, 600*cols), min(1000, 400*rows)]);
        
        % Distinct color palette for locations
        distinctColors = [
            0.00, 0.45, 0.74;  % Blue
            0.85, 0.33, 0.10;  % Red-Orange
            0.93, 0.69, 0.13;  % Yellow
            0.49, 0.18, 0.56;  % Purple
            0.47, 0.67, 0.19;  % Green
            0.30, 0.75, 0.93;  % Cyan
            0.64, 0.08, 0.18;  % Dark Red
            0.97, 0.51, 0.75;  % Pink
            0.00, 0.50, 0.00;  % Dark Green
            0.75, 0.75, 0.00;  % Olive
            0.00, 0.75, 0.75;  % Teal
            0.75, 0.00, 0.75;  % Magenta
            0.25, 0.25, 0.25;  % Dark Gray
            1.00, 0.50, 0.00;  % Orange
            0.00, 0.00, 0.55   % Navy Blue
        ];
        
        % Plot each patient in its own subplot
        for p = 1:numPatients
            subplot(rows, cols, p);
            hold on;
            
            % Get data for this patient
            dataRows = cleanData(3:end, :);
            
            % Filter by this patient
            if patientIDCol > 0
                keepRows = false(size(dataRows, 1), 1);
                for i = 1:size(dataRows, 1)
                    val = dataRows{i, patientIDCol};
                    numVal = toNumber(val);
                    if ~isnan(numVal) && numVal == patientNums(p)
                        keepRows(i) = true;
                    end
                end
                dataRows = dataRows(keepRows, :);
            end
            
            % Sort chronologically so plot reads left-to-right in time
            dataRows = sortRowsByTime(dataRows, headers);
            
            % Get unique catheter places for this patient
            catheterPlaces = {};
            for i = 1:size(dataRows, 1)
                val = dataRows{i, catheterPlaceCol};
                if ischar(val) && ~isempty(val)
                    catheterPlaces{end+1} = val;
                end
            end
            uniquePlaces = unique(catheterPlaces);
            
            % Assign colors to locations
            numPlaces = length(uniquePlaces);
            if numPlaces > 0
                colors = distinctColors(mod(0:numPlaces-1, size(distinctColors, 1)) + 1, :);
                
                % Collect all data: for each row, get its place, time label, and value
                % Use row index for x-positioning, then label x-axis with time labels
                allRowValues = [];     % parameter values
                allRowIndices = [];    % row indices (x positions)
                allRowColors = [];     % color per point
                allRowPlaces = {};     % place name per point
                
                % Build a place-to-color map
                placeColorMap = containers.Map();
                for loc = 1:numPlaces
                    placeColorMap(uniquePlaces{loc}) = loc;
                end
                
                for i = 1:size(dataRows, 1)
                    placeVal = dataRows{i, catheterPlaceCol};
                    if ischar(placeVal) || isstring(placeVal)
                        placeStr = char(placeVal);
                        if ~isempty(placeStr) && isKey(placeColorMap, placeStr)
                            paramVal = dataRows{i, colIdx};
                            numVal = toNumber(paramVal);
                            if ~isnan(numVal)
                                allRowValues(end+1) = numVal;
                                allRowIndices(end+1) = i;
                                colorIdx = placeColorMap(placeStr);
                                allRowColors(end+1, :) = colors(colorIdx, :);
                                allRowPlaces{end+1} = placeStr;
                            end
                        end
                    end
                end
                
                % Plot each location as a separate series for legend
                for loc = 1:numPlaces
                    place = uniquePlaces{loc};
                    mask = strcmp(allRowPlaces, place);
                    locIndices = allRowIndices(mask);
                    locValues = allRowValues(mask);
                    
                    if ~isempty(locValues)
                        plot(locIndices, locValues, '-o', ...
                            'LineWidth', 2, 'MarkerSize', 7, ...
                            'Color', colors(loc, :), ...
                            'MarkerFaceColor', colors(loc, :), ...
                            'DisplayName', place);
                    end
                end
                
                % Set x-axis tick labels using shared helper
                [tickPts, tickLbls] = buildTimeAxisLabels(dataRows, headers);
                if ~isempty(tickPts)
                    xticks(tickPts);
                    xticklabels(tickLbls);
                    xtickangle(45);
                end
                
                legend('Location', 'best', 'FontSize', 9);
            else
                text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center', 'FontSize', 12);
            end
            
            hold off;
            
            xlabel('Time', 'FontSize', 10, 'FontWeight', 'bold');
            ylabel(paramWithUnit, 'FontSize', 10, 'FontWeight', 'bold');
            title(sprintf('Patient %d', patientNums(p)), 'FontSize', 12, 'FontWeight', 'bold');
            grid on;
        end
        
        sgtitle(sprintf('%s Across Sample Locations', paramWithUnit), ...
            'FontSize', 14, 'FontWeight', 'bold');
        
        statusLabel.Text = sprintf('Plotted %s for %d patients', selectedParam, length(patientNums));
        statusLabel.FontColor = [0 0.5 0];
    end
    function plotParameterAcrossPlacesSingle(patientLabel, selectedParam)
        % Original single-patient location comparison logic
        selectedPatient = patientLabel;
        
        paramIdx = 0;
        for i = 1:length(paramNames)
            if strcmp(paramNames{i}, selectedParam)
                paramIdx = i;
                break;
            end
        end
        
        if paramIdx == 0
            statusLabel.Text = 'Parameter not found';
            statusLabel.FontColor = [1 0 0];
            return;
        end
        
        if paramIndices(paramIdx) == 0
            statusLabel.Text = 'Parameter not in data';
            statusLabel.FontColor = [1 0 0];
            return;
        end
        
        colIdx = paramIndices(paramIdx);
        
        dataRows = cleanData(3:end, :);
        keepAllPatients = strcmp(selectedPatient, 'All Patients');
        
        if ~keepAllPatients
            if patientIDCol > 0
                patientNum = str2double(selectedPatient);
                keepRows = false(size(dataRows, 1), 1);
                for i = 1:size(dataRows, 1)
                    val = dataRows{i, patientIDCol};
                    numVal = toNumber(val);
                    if ~isnan(numVal)
                        if numVal == patientNum
                            keepRows(i) = true;
                        end
                    end
                end
                dataRows = dataRows(keepRows, :);
            end
        end
        
        % Sort chronologically so plot reads left-to-right in time
        dataRows = sortRowsByTime(dataRows, headers);
        
        catheterPlaceCol = 0;
        for i = 1:length(headers)
            if ~isempty(headers{i})
                try
                    headerStr = string(headers{i});
                    if contains(headerStr, 'First Name', 'IgnoreCase', true)
                        catheterPlaceCol = i;
                        break;
                    end
                catch
                end
            end
        end
        
        timeCol = 0;
        for i = 1:length(headers)
            if ~isempty(headers{i})
                try
                    headerStr = string(headers{i});
                    if contains(headerStr, 'Last Name', 'IgnoreCase', true)
                        timeCol = i;
                        break;
                    end
                catch
                end
            end
        end
        
        % Find actual Time column as fallback
        actualTimeCol = 0;
        for i = 1:length(headers)
            if ~isempty(headers{i})
                try
                    headerStr = string(headers{i});
                    if strcmpi(strtrim(headerStr), 'Time')
                        actualTimeCol = i;
                        break;
                    end
                catch
                end
            end
        end
        
        if catheterPlaceCol == 0
            statusLabel.Text = 'First Name column not found';
            statusLabel.FontColor = [1 0 0];
            return;
        end
        
        if timeCol == 0 && actualTimeCol == 0
            statusLabel.Text = 'No time column found (Last Name or Time)';
            statusLabel.FontColor = [1 0 0];
            return;
        end
        
        catheterPlaces = {};
        for i = 1:size(dataRows, 1)
            val = dataRows{i, catheterPlaceCol};
            if ischar(val)
                if ~isempty(val)
                    catheterPlaces{end+1} = val;
                end
            end
        end
        uniquePlaces = unique(catheterPlaces);
        
        % Create new figure
        placeFig = figure('Name', sprintf('%s Across Catheter Places - Patient %s', ...
            selectedParam, selectedPatient), ...
            'Position', [100 100 1400 800]);
        
        % Define highly distinct color palette for better visualization
        distinctColors = [
            0.00, 0.45, 0.74;  % Blue
            0.85, 0.33, 0.10;  % Red-Orange
            0.93, 0.69, 0.13;  % Yellow
            0.49, 0.18, 0.56;  % Purple
            0.47, 0.67, 0.19;  % Green
            0.30, 0.75, 0.93;  % Cyan
            0.64, 0.08, 0.18;  % Dark Red
            0.97, 0.51, 0.75;  % Pink
            0.00, 0.50, 0.00;  % Dark Green
            0.75, 0.75, 0.00;  % Olive
            0.00, 0.75, 0.75;  % Teal
            0.75, 0.00, 0.75;  % Magenta
            0.25, 0.25, 0.25;  % Dark Gray
            1.00, 0.50, 0.00;  % Orange
            0.00, 0.00, 0.55   % Navy Blue
        ];
        
        % Use distinct colors, repeat if more locations than colors
        numPlaces = length(uniquePlaces);
        colors = distinctColors(mod(0:numPlaces-1, size(distinctColors, 1)) + 1, :);
        
        hold on;
        
        for p = 1:length(uniquePlaces)
            place = uniquePlaces{p};
            
            % Filter data for this place
            placeValues = [];
            placeTimeLabels = {};
            placeIndices = [];
            
            for i = 1:size(dataRows, 1)
                placeVal = dataRows{i, catheterPlaceCol};
                if ischar(placeVal) || isstring(placeVal)
                    if strcmp(char(placeVal), place)
                        paramVal = dataRows{i, colIdx};
                        numVal = toNumber(paramVal);
                        if ~isnan(numVal)
                            if numVal ~= 0
                                placeValues(end+1) = numVal;
                                placeIndices(end+1) = i;
                                
                                timeVal = dataRows{i, timeCol};
                                if ischar(timeVal) || isstring(timeVal)
                                    if ~isempty(timeVal)
                                        placeTimeLabels{end+1} = char(timeVal);
                                    else
                                        placeTimeLabels{end+1} = sprintf('%d', i);
                                    end
                                else
                                    placeTimeLabels{end+1} = sprintf('%d', i);
                                end
                            end
                        end
                    end
                end
            end
            
            if ~isempty(placeValues)
                plot(placeIndices, placeValues, '-o', ...
                    'LineWidth', 2.5, 'MarkerSize', 8, ...
                    'Color', colors(p, :), ...
                    'DisplayName', place);
            end
        end
        
        hold off;
        
        % Set x-axis tick labels using shared helper
        [tickPts, tickLbls] = buildTimeAxisLabels(dataRows, headers);
        if ~isempty(tickPts)
            xticks(tickPts);
            xticklabels(tickLbls);
            xtickangle(45);
        end
        
        xlabel('Time', 'FontSize', 12, 'FontWeight', 'bold');
        paramWithUnit = getParamWithUnit(selectedParam);
        ylabel(paramWithUnit, 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('%s Across Catheter Places Over Time - Patient %s', ...
            paramWithUnit, selectedPatient), ...
            'FontSize', 14, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 11);
        grid on;
        
        statusLabel.Text = sprintf('Plotted %s across %d places', selectedParam, length(uniquePlaces));
        statusLabel.FontColor = [0 0.5 0];
    end  % end plotParameterAcrossPlacesSingle
    
    function calculateStatistics()
        if isempty(cleanData)
            return;
        end
        
        % Create statistics dialog with better sizing
        statsFig = uifigure('Name', 'Clinical Statistical Analysis', 'Position', [100 100 1000 750]);
        
        statsGrid = uigridlayout(statsFig, [6 2]);
        statsGrid.RowHeight = {'fit', 180, 'fit', 95, 'fit', '1x'};
        statsGrid.ColumnWidth = {'1x', '1x'};
        statsGrid.Padding = [15 15 15 15];
        statsGrid.RowSpacing = 10;
        statsGrid.ColumnSpacing = 15;
        
        titleLabel = uilabel(statsGrid, 'Text', 'Clinical Statistical Analysis', ...
            'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'FontColor', [0.2 0.4 0.7]);
        titleLabel.Layout.Row = 1;
        titleLabel.Layout.Column = [1 2];
        
        patientPanel = uipanel(statsGrid, 'Title', 'Select Patients');
        patientPanel.Layout.Row = 2;
        patientPanel.Layout.Column = 1;
        
        patientPanelGrid = uigridlayout(patientPanel, [3 1]);
        patientPanelGrid.RowHeight = {'fit', '1x', 'fit'};
        patientPanelGrid.Padding = [5 5 5 5];
        
        instructionLabel = uilabel(patientPanelGrid, ...
            'Text', 'Ctrl+Click to select multiple from list:', ...
            'FontSize', 9, 'FontColor', [0.4 0.4 0.4]);
        
        if patientIDCol > 0
            dataRows = cleanData(3:end, :);
            patientIDs = [];
            for i = 1:size(dataRows, 1)
                val = dataRows{i, patientIDCol};
                numVal = toNumber(val);
                if ~isnan(numVal)
                    patientIDs(end+1) = numVal;
                end
            end
            uniquePatients = unique(patientIDs);
            patientList = arrayfun(@num2str, uniquePatients, 'UniformOutput', false);
        else
            patientList = {'All'};
        end
        
        patientsListBox = uilistbox(patientPanelGrid, ...
            'Items', patientList, ...
            'Multiselect', 'on', ...
            'ValueChangedFcn', @(src, event) updateTimePointsForPatients());
        
        manualEntryField = uieditfield(patientPanelGrid, 'text', ...
            'Placeholder', 'Or type IDs here (comma-separated): 2025050,2025051', ...
            'FontSize', 9, ...
            'ValueChangedFcn', @(src, event) updateTimePointsForPatients());
        
        rightPanel = uipanel(statsGrid, 'Title', 'Filters');
        rightPanel.Layout.Row = 2;
        rightPanel.Layout.Column = 2;
        rightPanelGrid = uigridlayout(rightPanel, [4 1]);
        rightPanelGrid.RowHeight = {'fit', 'fit', 'fit', 'fit'};
        rightPanelGrid.Padding = [10 10 10 10];
        
        uilabel(rightPanelGrid, 'Text', 'Catheter Place:', 'FontWeight', 'bold');
        statsCatheterDropdown = uidropdown(rightPanelGrid, 'Items', catheterPlaceDropdown.Items);
        
        uilabel(rightPanelGrid, 'Text', 'Time Point:', 'FontWeight', 'bold');
        
        timeCol = 0;
        for i = 1:length(headers)
            if ~isempty(headers{i})
                try
                    headerStr = string(headers{i});
                    if contains(headerStr, 'Last Name', 'IgnoreCase', true)
                        timeCol = i;
                        break;
                    end
                catch
                end
            end
        end
        
        timePoints = {'All'};
        if timeCol > 0
            dataRows = cleanData(3:end, :);
            for i = 1:size(dataRows, 1)
                timeVal = dataRows{i, timeCol};
                if ischar(timeVal) || isstring(timeVal)
                    if ~isempty(timeVal)
                        timePoints{end+1} = char(timeVal);
                    end
                end
            end
            timePoints = unique(timePoints);
        end
        
        statsTimeDropdown = uidropdown(rightPanelGrid, 'Items', timePoints);
        
        function updateTimePointsForPatients()
            selectedPatientIDs = {};
            
            manualEntry = strtrim(manualEntryField.Value);
            if ~isempty(manualEntry)
                manualIDs = strsplit(manualEntry, ',');
                for i = 1:length(manualIDs)
                    id = strtrim(manualIDs{i});
                    if ~isempty(id)
                        selectedPatientIDs{end+1} = id;
                    end
                end
            else
                selectedPatientIDs = patientsListBox.Value;
            end
            
            if isempty(selectedPatientIDs)
                statsTimeDropdown.Items = {'All'};
                statsCatheterDropdown.Items = {'All Places'};
                return;
            end
            
            updateCatheterPlacesForPatients(selectedPatientIDs);
            
            updateTimePoints(selectedPatientIDs);
        end
        
        function updateCatheterPlacesForPatients(selectedPatientIDs)
            catheterPlaceCol = 0;
            for i = 1:length(headers)
                if ~isempty(headers{i})
                    try
                        headerStr = string(headers{i});
                        if contains(headerStr, 'First Name', 'IgnoreCase', true)
                            catheterPlaceCol = i;
                            break;
                        end
                    catch
                    end
                end
            end
            
            if catheterPlaceCol > 0 && patientIDCol > 0
                dataRows = cleanData(3:end, :);
                
                patientCatheterSets = {};
                for p = 1:length(selectedPatientIDs)
                    patientID = str2double(selectedPatientIDs{p});
                    patientCatheters = {};
                    
                    for i = 1:size(dataRows, 1)
                        patVal = dataRows{i, patientIDCol};
                        numVal = toNumber(patVal);
                        
                        if ~isnan(numVal) && numVal == patientID
                            cathVal = dataRows{i, catheterPlaceCol};
                            if ischar(cathVal) || isstring(cathVal)
                                cathStr = strtrim(char(cathVal));
                                if ~isempty(cathStr)
                                    patientCatheters{end+1} = cathStr;
                                end
                            end
                        end
                    end
                    
                    if ~isempty(patientCatheters)
                        patientCatheterSets{end+1} = unique(patientCatheters);
                    end
                end
                
                if ~isempty(patientCatheterSets)
                    commonCatheters = patientCatheterSets{1};
                    
                    for p = 2:length(patientCatheterSets)
                        commonCatheters = intersect(commonCatheters, patientCatheterSets{p});
                    end
                    
                    if ~isempty(commonCatheters)
                        statsCatheterDropdown.Items = ['All Places', commonCatheters];
                    else
                        statsCatheterDropdown.Items = {'All Places', '(No common catheter places)'};
                    end
                else
                    statsCatheterDropdown.Items = {'All Places'};
                end
                
                statsCatheterDropdown.Value = 'All Places';
            end
        end
        
        function updateTimePoints(selectedPatientIDs)
            if timeCol > 0 && patientIDCol > 0
                dataRows = cleanData(3:end, :);
                
                patientTimeSets = {};
                for p = 1:length(selectedPatientIDs)
                    patientID = str2double(selectedPatientIDs{p});
                    patientTimes = {};
                    
                    for i = 1:size(dataRows, 1)
                        patVal = dataRows{i, patientIDCol};
                        numVal = toNumber(patVal);
                        
                        if ~isnan(numVal) && numVal == patientID
                            timeVal = dataRows{i, timeCol};
                            if ischar(timeVal) || isstring(timeVal)
                                timeStr = char(timeVal);
                                if ~isempty(timeStr)
                                    patientTimes{end+1} = timeStr;
                                end
                            end
                        end
                    end
                    
                    if ~isempty(patientTimes)
                        patientTimeSets{end+1} = unique(patientTimes);
                    end
                end
                
                if ~isempty(patientTimeSets)
                    commonTimes = patientTimeSets{1};
                    
                    for p = 2:length(patientTimeSets)
                        commonTimes = intersect(commonTimes, patientTimeSets{p});
                    end
                    
                    if ~isempty(commonTimes)
                        statsTimeDropdown.Items = ['All', commonTimes];
                    else
                        statsTimeDropdown.Items = {'All', '(No common time points)'};
                    end
                else
                    statsTimeDropdown.Items = {'All'};
                end
                
                statsTimeDropdown.Value = 'All';
            end
        end
        
        paramLabel = uilabel(statsGrid, 'Text', 'Select Parameters:', ...
            'FontWeight', 'bold', 'FontSize', 12);
        paramLabel.Layout.Row = 3;
        paramLabel.Layout.Column = [1 2];
        
        paramPanel = uipanel(statsGrid);
        paramPanel.Layout.Row = 4;
        paramPanel.Layout.Column = [1 2];
        paramPanelGrid = uigridlayout(paramPanel, [ceil(length(paramNames)/4), 4]);
        paramPanelGrid.RowHeight = repmat({'fit'}, 1, ceil(length(paramNames)/4));
        paramPanelGrid.ColumnWidth = {'1x', '1x', '1x', '1x'};
        
        paramCheckboxes = cell(1, length(paramNames));
        for i = 1:length(paramNames)
            paramCheckboxes{i} = uicheckbox(paramPanelGrid, 'Text', paramNames{i}, 'Value', false);
        end
        
        % Calculate button - spans both columns
        calculateBtn = uibutton(statsGrid, 'Text', '📊 Calculate Statistics', ...
            'ButtonPushedFcn', @(btn,event) performCalculation(), ...
            'BackgroundColor', [0.2 0.8 0.4], 'FontColor', 'white', ...
            'FontSize', 14, 'FontWeight', 'bold');
        calculateBtn.Layout.Row = 5;
        calculateBtn.Layout.Column = [1 2];
        
        resultsArea = uitextarea(statsGrid, 'Editable', 'off', ...
            'Value', {'Select patients, filters, and parameters, then click Calculate Statistics'}, ...
            'FontSize', 11, 'FontName', 'Courier New');
        resultsArea.Layout.Row = 6;
        resultsArea.Layout.Column = [1 2];
        
        function performCalculation()
            selectedPatientsList = {};
            
            manualEntry = strtrim(manualEntryField.Value);
            if ~isempty(manualEntry)
                manualIDs = strsplit(manualEntry, ',');
                for i = 1:length(manualIDs)
                    id = strtrim(manualIDs{i});
                    if ~isempty(id)
                        selectedPatientsList{end+1} = id;
                    end
                end
            else
                selectedPatientsList = patientsListBox.Value;
            end
            
            if isempty(selectedPatientsList)
                resultsArea.Value = {'Please select patients from the list or type IDs in the text field'};
                return;
            end
            
            selectedPlace = statsCatheterDropdown.Value;
            
            selectedTime = statsTimeDropdown.Value;
            
            selectedParams = {};
            selectedParamIndices = [];
            for i = 1:length(paramCheckboxes)
                if paramCheckboxes{i}.Value
                    if paramIndices(i) > 0
                        selectedParams{end+1} = paramNames{i};
                        selectedParamIndices(end+1) = paramIndices(i);
                    end
                end
            end
            
            if isempty(selectedParams)
                resultsArea.Value = {'Please select at least one parameter'};
                return;
            end
            
            % Filter data
            dataRows = cleanData(3:end, :);
            
            % Filter by patients
            keepRows = false(size(dataRows, 1), 1);
            for i = 1:size(dataRows, 1)
                patVal = dataRows{i, patientIDCol};
                numVal = toNumber(patVal);
                if ~isnan(numVal)
                    patStr = num2str(numVal);
                    for j = 1:length(selectedPatientsList)
                        selectedPatID = selectedPatientsList{j};
                        if strcmp(patStr, selectedPatID)
                            keepRows(i) = true;
                            break;
                        end
                    end
                end
            end
            dataRows = dataRows(keepRows, :);
            
            % Filter by catheter place
            if ~strcmp(selectedPlace, 'All Places')
                catheterPlaceCol = 0;
                for i = 1:length(headers)
                    if ~isempty(headers{i})
                        try
                            headerStr = string(headers{i});
                            if contains(headerStr, 'First Name', 'IgnoreCase', true)
                                catheterPlaceCol = i;
                                break;
                            end
                        catch
                        end
                    end
                end
                
                if catheterPlaceCol > 0
                    keepRows = false(size(dataRows, 1), 1);
                    for i = 1:size(dataRows, 1)
                        val = dataRows{i, catheterPlaceCol};
                        if ischar(val)
                            if strcmp(val, selectedPlace)
                                keepRows(i) = true;
                            end
                        end
                    end
                    dataRows = dataRows(keepRows, :);
                end
            end
            
            % Filter by time
            if ~strcmp(selectedTime, 'All')
                if timeCol > 0
                    keepRows = false(size(dataRows, 1), 1);
                    for i = 1:size(dataRows, 1)
                        val = dataRows{i, timeCol};
                        if ischar(val)
                            if strcmp(val, selectedTime)
                                keepRows(i) = true;
                            end
                        end
                    end
                    dataRows = dataRows(keepRows, :);
                end
            end
            
            if size(dataRows, 1) == 0
                resultsArea.Value = {'No data found with the selected criteria'};
                return;
            end
            
            % Calculate statistics
            results = {};
            results{end+1} = '========================================';
            results{end+1} = 'STATISTICAL ANALYSIS RESULTS';
            results{end+1} = '========================================';
            results{end+1} = '';
            results{end+1} = sprintf('Patients: %s', strjoin(selectedPatientsList, ', '));
            results{end+1} = sprintf('Catheter Place: %s', selectedPlace);
            results{end+1} = sprintf('Time Point: %s', selectedTime);
            results{end+1} = sprintf('Sample Size: %d measurements', size(dataRows, 1));
            results{end+1} = '';
            results{end+1} = '========================================';
            
            for p = 1:length(selectedParams)
                colIdx = selectedParamIndices(p);
                
                values = [];
                for i = 1:size(dataRows, 1)
                    val = dataRows{i, colIdx};
                    numVal = toNumber(val);
                    if ~isnan(numVal)
                        values(end+1) = numVal;
                    end
                end
                
                if ~isempty(values)
                    meanVal = mean(values);
                    stdVal = std(values);
                    minVal = min(values);
                    maxVal = max(values);
                    medianVal = median(values);
                    
                    results{end+1} = '';
                    results{end+1} = sprintf('--- %s ---', selectedParams{p});
                    results{end+1} = sprintf('  Mean:   %.2f', meanVal);
                    results{end+1} = sprintf('  SD:     %.2f', stdVal);
                    results{end+1} = sprintf('  Median: %.2f', medianVal);
                    results{end+1} = sprintf('  Min:    %.2f', minVal);
                    results{end+1} = sprintf('  Max:    %.2f', maxVal);
                    results{end+1} = sprintf('  N:      %d', length(values));
                else
                    results{end+1} = '';
                    results{end+1} = sprintf('--- %s ---', selectedParams{p});
                    results{end+1} = '  No valid data';
                end
            end
            
            results{end+1} = '';
            results{end+1} = '========================================';
            
            resultsArea.Value = results;
            
            createStatsChart();
        end
        
        function createStatsChart()
            selectedPatientsList = {};
            
            manualEntry = strtrim(manualEntryField.Value);
            if ~isempty(manualEntry)
                manualIDs = strsplit(manualEntry, ',');
                for i = 1:length(manualIDs)
                    id = strtrim(manualIDs{i});
                    if ~isempty(id)
                        selectedPatientsList{end+1} = id;
                    end
                end
            else
                selectedPatientsList = patientsListBox.Value;
            end
            
            % Create figure with boxplots
            chartFig = figure('Name', 'Statistics Boxplot', 'Position', [150 150 1200 700]);
            
            selectedParams = {};
            allValues = {};
            
            dataRows = cleanData(3:end, :);
            
            keepRows = false(size(dataRows, 1), 1);
            for i = 1:size(dataRows, 1)
                patVal = dataRows{i, patientIDCol};
                numVal = toNumber(patVal);
                if ~isnan(numVal)
                    patStr = num2str(numVal);
                    for j = 1:length(selectedPatientsList)
                        selectedPatID = selectedPatientsList{j};
                        if strcmp(patStr, selectedPatID)
                            keepRows(i) = true;
                            break;
                        end
                    end
                end
            end
            dataRows = dataRows(keepRows, :);
            
            if ~strcmp(statsCatheterDropdown.Value, 'All Places')
                catheterPlaceCol = 0;
                for i = 1:length(headers)
                    if ~isempty(headers{i})
                        try
                            if contains(string(headers{i}), 'First Name', 'IgnoreCase', true)
                                catheterPlaceCol = i;
                                break;
                            end
                        catch
                        end
                    end
                end
                
                if catheterPlaceCol > 0
                    keepRows = false(size(dataRows, 1), 1);
                    for i = 1:size(dataRows, 1)
                        val = dataRows{i, catheterPlaceCol};
                        if ischar(val)
                            if strcmp(val, statsCatheterDropdown.Value)
                                keepRows(i) = true;
                            end
                        end
                    end
                    dataRows = dataRows(keepRows, :);
                end
            end
            
            if ~strcmp(statsTimeDropdown.Value, 'All')
                if timeCol > 0
                    keepRows = false(size(dataRows, 1), 1);
                    for i = 1:size(dataRows, 1)
                        val = dataRows{i, timeCol};
                        if ischar(val)
                            if strcmp(val, statsTimeDropdown.Value)
                                keepRows(i) = true;
                            end
                        end
                    end
                    dataRows = dataRows(keepRows, :);
                end
            end
            
            for i = 1:length(paramCheckboxes)
                if paramCheckboxes{i}.Value
                    if paramIndices(i) > 0
                        colIdx = paramIndices(i);
                        
                        values = [];
                        for j = 1:size(dataRows, 1)
                            val = dataRows{j, colIdx};
                            numVal = toNumber(val);
                            if ~isnan(numVal)
                                values(end+1) = numVal;
                            end
                        end
                        
                        if ~isempty(values)
                            selectedParams{end+1} = paramNames{i};
                            allValues{end+1} = values;
                        end
                    end
                end
            end
            
            if ~isempty(allValues)
                maxLen = max(cellfun(@length, allValues));
                boxData = nan(maxLen, length(allValues));
                
                for i = 1:length(allValues)
                    vals = allValues{i};
                    boxData(1:length(vals), i) = vals;
                end
                
                % Create boxplot with proper group labels
                groups = [];
                allData = [];
                groupLabels = {};
                
                for i = 1:length(allValues)
                    vals = allValues{i};
                    allData = [allData; vals(:)];
                    groups = [groups; repmat(i, length(vals), 1)];
                    groupLabels{i} = selectedParams{i};
                end
                
                % Create boxplot
                boxplot(allData, groups, 'Labels', groupLabels, 'Colors', 'b', ...
                    'Widths', 0.5, 'Symbol', 'o');
                
                xtickangle(45);
                ylabel('Values', 'FontSize', 14, 'FontWeight', 'bold');
                title(sprintf('Boxplot - Patients: %s, Place: %s, Time: %s', ...
                    strjoin(selectedPatientsList, ','), ...
                    statsCatheterDropdown.Value, statsTimeDropdown.Value), ...
                    'FontSize', 14, 'FontWeight', 'bold');
                grid on;
                
                hold on;
                means = cellfun(@mean, allValues);
                plot(1:length(means), means, 'rd', 'MarkerSize', 10, ...
                    'MarkerFaceColor', 'r', 'DisplayName', 'Mean');
                legend('show', 'Location', 'best');
                hold off;
            end
        end
    end

end
