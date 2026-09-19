function launch_PatLogGUI(varargin)
%LAUNCH_PATLOGGUI Entry point for the compiled standalone application.
%
%   A compiled application exits as soon as its entry function returns, so
%   this blocks until the user closes the app window. Any command-line
%   argument is passed through as the PatLog export to open on startup.

    PatLogGUI(varargin{:});

    figs = findall(0, 'Type', 'figure');
    keep = arrayfun(@(f) contains(string(f.Name), 'PatLog Blood Gas Analyzer'), figs);
    if any(keep)
        waitfor(figs(find(keep, 1)));
    end
end
