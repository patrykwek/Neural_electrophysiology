%% ========= SESSION PROGRESSION PLOT: MEAN CUE DECODING IN SIGNIFICANT INTERVAL =========
% Requires these variables from your decoding script:
%   AccMat        [nBins x nSessions]
%   xPlot         [1 x nBins]
%   obsClusters   cell array of cluster bin indices
%   obsClusterP   p-value per observed cluster
%   nSessions
%   outDir
%
% This plot answers:
%   Does cue decoding strength gradually increase across sessions?

fprintf('\n=== SESSION PROGRESSION PLOT (SIGNIFICANT INTERVAL) ===\n');

% ---- collect significant clusters from cluster permutation result ----
sigClusterIdx = find(obsClusterP < 0.05);

if isempty(sigClusterIdx)
    warning('No significant clusters found (obsClusterP < 0.05). Session progression plot not created.');
else
    % If multiple significant clusters exist, use all of them together
    sigBinsMask = false(1, size(AccMat,1));
    for iC = 1:numel(sigClusterIdx)
        sigBinsMask(obsClusters{sigClusterIdx(iC)}) = true;
    end

    sigBins = find(sigBinsMask);
    sigTimeStart = xPlot(sigBins(1));
    sigTimeEnd   = xPlot(sigBins(end));

    fprintf('Using significant interval from %.3f to %.3f s\n', sigTimeStart, sigTimeEnd);

    % ---- mean decoding within significant interval for each session ----
    sessSigMeanAcc = nan(1, nSessions);
    for sIdx = 1:nSessions
        sessSigMeanAcc(sIdx) = mean(AccMat(sigBins, sIdx), 'omitnan');
    end

    sessionNums = 1:nSessions;
    validSess = isfinite(sessSigMeanAcc);

    if nnz(validSess) < 2
        warning('Not enough valid sessions to plot session progression.');
    else
        % ---- linear fit ----
        xFit = sessionNums(validSess)';
        yFit = sessSigMeanAcc(validSess)';

        pLin = polyfit(xFit, yFit, 1);
        yHat = polyval(pLin, xFit);

        SSres = sum((yFit - yHat).^2);
        SStot = sum((yFit - mean(yFit)).^2);
        R2 = 1 - SSres / SStot;

        % Pearson correlation
        [rVal, pVal] = corr(xFit, yFit, 'Type', 'Pearson');

        % ---- figure style matched to your other plots ----
        figProg = figure('Color','w','Position',[120,120,700,500]);
        axProg = axes(figProg); hold(axProg,'on');

        plot(axProg, sessionNums(validSess), sessSigMeanAcc(validSess), 'o', ...
            'MarkerSize', 10, ...
            'MarkerFaceColor', greyDot, ...
            'MarkerEdgeColor', greyDot, ...
            'LineWidth', scatterEdgeLW);

        plot(axProg, xFit, yHat, '-', ...
            'Color', [0 0 0], ...
            'LineWidth', traceLW);

        xlabel(axProg, 'Session', 'FontSize', labelFontSize);
        ylabel(axProg, 'Mean decoding accuracy', 'FontSize', labelFontSize);
        title(axProg, sprintf('Session progression of cue decoding (%.2f to %.2f s)', sigTimeStart, sigTimeEnd), ...
            'FontWeight', 'bold', 'FontSize', titleFontSize);

        set(axProg, 'FontSize', tickFontSize, ...
            'TickDir','out', ...
            'LineWidth', axesTickLW, ...
            'TickLength', axesTickLen, ...
            'Box','off', ...
            'Layer','top');

        xlim(axProg, [1 nSessions]);

        % optional chance line for 3-class decoding
        yline(axProg, 1/3, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 2);

        % text annotation
        xl = xlim(axProg);
        yl = ylim(axProg);
        xTxt = xl(1) + 0.03 * (xl(2)-xl(1));
        yTxt = yl(2) - 0.08 * (yl(2)-yl(1));

        txtStr = sprintf('r = %.3f, p = %.4g\nR^2 = %.3f\nslope = %.5f', rVal, pVal, R2, pLin(1));
        text(axProg, xTxt, yTxt, txtStr, ...
            'FontSize', 18, ...
            'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'left');

        % ---- save figure ----
        outProgBase = fullfile(outDir, ...
            sprintf('Sessions01_%02d_CueTypeDecoding_SessionProgression_SigInterval_%0.3fto%0.3fs', ...
            nSessions, sigTimeStart, sigTimeEnd));

        set(figProg, 'PaperPositionMode', 'auto');
        exportgraphics(figProg, [outProgBase '.png'], 'Resolution', 300, 'BackgroundColor', 'white');
        savefig(figProg, [outProgBase '.fig']);
        try
            saveas(figProg, [outProgBase '.svg']);
        catch
            print(figProg, [outProgBase '.svg'], '-dsvg');
        end

        fprintf('Saved:\n  %s.png\n  %s.fig\n  %s.svg\n', outProgBase, outProgBase, outProgBase);

        % ---- save underlying values ----
        outProgMat = [outProgBase '.mat'];
        save(outProgMat, ...
            'sessSigMeanAcc', ...
            'sessionNums', ...
            'validSess', ...
            'sigBins', ...
            'sigBinsMask', ...
            'sigTimeStart', ...
            'sigTimeEnd', ...
            'xFit', ...
            'yFit', ...
            'yHat', ...
            'pLin', ...
            'R2', ...
            'rVal', ...
            'pVal', ...
            '-v7.3');

        fprintf('Saved session progression data:\n  %s\n', outProgMat);
    end
end