function launch_ABL_CDI_Analyzer()
%LAUNCH_ABL_CDI_ANALYZER Entry point for the compiled standalone application.
%
%   A compiled application exits as soon as its entry function returns, so
%   this blocks until the user closes the app window. ABL_CDI_Analyzer is a
%   classdef, which cannot itself be a compiler entry point.

    app = ABL_CDI_Analyzer;
    waitfor(app.UIFigure);
end
