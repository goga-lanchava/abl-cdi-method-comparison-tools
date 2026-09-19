function make_branding(outDir)
%MAKE_BRANDING Generate the icons and splash screens used by the installers.
%
%   Two related tools, so one visual language: a deep clinical teal ground
%   with a light motif that states what each tool does.
%
%     ABL-CDI Analyzer - a continuous CDI trace with discrete ABL reference
%                        draws sitting on it, which is the whole premise of
%                        the method comparison.
%     PatLogGUI        - a data grid with a trend lifted out of it, which is
%                        the import-clean-explore path.
%
%   make_branding(outDir) writes the PNGs into outDir (default: this folder).

    if nargin < 1 || isempty(outDir)
        outDir = fileparts(mfilename('fullpath'));
    end
    if ~exist(outDir, 'dir'); mkdir(outDir); end

    TEAL_DARK = [0.031 0.271 0.318];   % #082245-ish deep teal
    TEAL_MID  = [0.051 0.420 0.471];   % #0d6b78
    AMBER     = [0.953 0.694 0.298];   % discrete reference draws
    PAPER     = [0.945 0.969 0.976];

    drawIcon(fullfile(outDir,'abl_icon.png'),    @motifAnalyzer, 512, TEAL_DARK, TEAL_MID, AMBER, PAPER);
    drawIcon(fullfile(outDir,'patlog_icon.png'), @motifPatLog,   512, TEAL_DARK, TEAL_MID, AMBER, PAPER);

    drawSplash(fullfile(outDir,'abl_splash.png'), 'ABL-CDI Analyzer', ...
        'Method comparison and correction', @motifAnalyzer, TEAL_DARK, TEAL_MID, AMBER, PAPER);
    drawSplash(fullfile(outDir,'patlog_splash.png'), 'PatLogGUI', ...
        'Blood gas data exploration', @motifPatLog, TEAL_DARK, TEAL_MID, AMBER, PAPER);

    fprintf('branding written to %s\n', outDir);
end

% ------------------------------------------------------------------------

function drawIcon(file, motifFcn, px, cDark, cMid, cAmber, cPaper)
    f = figure('Visible','off','Units','pixels','Position',[60 60 px px],'Color','w');
    ax = axes(f,'Position',[0 0 1 1]); hold(ax,'on');
    paintGround(ax, cDark, cMid);
    motifFcn(ax, cPaper, cAmber, 1.0, [0.08 0.92], [0.18 0.82]);
    axis(ax,[0 1 0 1]); axis(ax,'off');

    A = renderToImage(f, [px px]);
    imwrite(A, file, 'Alpha', roundedAlpha(px, px, round(px*0.18)));
    close(f);
end

function drawSplash(file, titleTxt, subTxt, motifFcn, cDark, cMid, cAmber, cPaper)
    W = 800; H = 500;
    f = figure('Visible','off','Units','pixels','Position',[60 60 W H],'Color','w');
    ax = axes(f,'Position',[0 0 1 1]); hold(ax,'on');
    paintGround(ax, cDark, cMid);

    % motif occupies the right panel only, so it never crosses the wordmark
    motifFcn(ax, cPaper, cAmber, 0.85, [0.50 0.94], [0.24 0.76]);

    text(ax, 0.07, 0.66, titleTxt, 'Color', cPaper, 'FontSize', 30, ...
        'FontWeight','bold', 'FontName','Helvetica', 'Units','normalized');
    text(ax, 0.072, 0.54, subTxt, 'Color', cPaper, 'FontSize', 14, ...
        'FontName','Helvetica', 'Units','normalized');
    plot(ax, [0.07 0.26], [0.475 0.475], '-', 'Color', cAmber, 'LineWidth', 3);
    text(ax, 0.072, 0.36, 'v1.1.0', 'Color', cPaper, 'FontSize', 11, ...
        'FontName','Helvetica', 'Units','normalized');
    text(ax, 0.072, 0.27, 'Czech Technical University in Prague', ...
        'Color', cPaper, 'FontSize', 10, 'FontName','Helvetica', 'Units','normalized');

    axis(ax,[0 1 0 1]); axis(ax,'off');
    imwrite(renderToImage(f, [H W]), file);
    close(f);
end

% ------------------------------------------------------------------------

function paintGround(ax, cDark, cMid)
    % vertical gradient, dark at the bottom
    ny = 256;
    g = zeros(ny, 2, 3);
    for k = 1:3
        g(:,:,k) = repmat(linspace(cDark(k), cMid(k), ny)', 1, 2);
    end
    image(ax, 'XData',[0 1], 'YData',[0 1], 'CData', flipud(g));
end

function motifAnalyzer(ax, cPaper, cAmber, alpha, xr, yr)
    % continuous CDI trace with discrete ABL draws on it
    w = diff(xr); h = diff(yr); yc = mean(yr);
    u = linspace(0, 1, 300);
    shape = @(v) yc + h*0.26*sin(2*pi*v/0.65) + h*0.09*sin(2*pi*v/0.21);
    plot(ax, xr(1) + w*u, shape(u), '-', 'Color', [cPaper alpha], 'LineWidth', 6);
    plot(ax, xr, [yc yc], ':', 'Color', [cPaper alpha*0.45], 'LineWidth', 2.5);
    us = linspace(0.08, 0.92, 5);
    scatter(ax, xr(1) + w*us, shape(us) + h*0.07, 520*w*2.2, 'filled', ...
        'MarkerFaceColor', cAmber, 'MarkerFaceAlpha', alpha, 'MarkerEdgeColor','none');
end

function motifPatLog(ax, cPaper, cAmber, alpha, xr, yr)
    % data grid with a trend lifted out of it
    w = diff(xr); h = diff(yr);
    for gx = linspace(xr(1)+w*0.08, xr(2)-w*0.08, 5)
        plot(ax, [gx gx], yr, '-', 'Color', [cPaper alpha*0.30], 'LineWidth', 2);
    end
    for gy = linspace(yr(1)+h*0.08, yr(2)-h*0.08, 4)
        plot(ax, xr, [gy gy], '-', 'Color', [cPaper alpha*0.30], 'LineWidth', 2);
    end
    us = linspace(0.08, 0.92, 5);
    vs = [0.18 0.42 0.33 0.72 0.90];
    plot(ax, xr(1) + w*us, yr(1) + h*vs, '-', 'Color', [cPaper alpha], 'LineWidth', 6);
    scatter(ax, xr(1) + w*us, yr(1) + h*vs, 420*w*2.2, 'filled', ...
        'MarkerFaceColor', cAmber, 'MarkerFaceAlpha', alpha, 'MarkerEdgeColor','none');
end

% ------------------------------------------------------------------------

function A = renderToImage(f, targetSize)
    tmp = [tempname '.png'];
    exportgraphics(f, tmp, 'Resolution', 150, 'BackgroundColor','current');
    A = imread(tmp);
    delete(tmp);
    A = imresize(A, targetSize);
end

function alphaMask = roundedAlpha(h, w, r)
    alphaMask = ones(h, w);
    [X, Y] = meshgrid(1:w, 1:h);
    corners = [r+1 r+1; w-r h-r; r+1 h-r; w-r r+1];
    inCorner = (X < r+1 & Y < r+1) | (X > w-r & Y > h-r) | ...
               (X < r+1 & Y > h-r) | (X > w-r & Y < r+1);
    d = inf(h, w);
    for c = 1:4
        d = min(d, hypot(X - corners(c,1), Y - corners(c,2)));
    end
    alphaMask(inCorner & d > r) = 0;
    soft = inCorner & d > r-1.5 & d <= r;
    alphaMask(soft) = max(0, min(1, r - d(soft)));
end
