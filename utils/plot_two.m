function fig = plot_two(x, y1, y2, lgnd, options)
%PLOT_TWO Summary of this function goes here
%   Detailed explanation goes here
    arguments
        x
        y1
        y2
        lgnd.ColorLegend
        lgnd.LineLegend
        lgnd.ColorLocation
        lgnd.LineLocation
        options.LineWidth
        options.Color
        options.XLabel
        options.YLabel
        options.Title
        options.XLim (1,:) double = NaN
        options.YLim (1,:) double = NaN
        options.Position (1,:) double = NaN
    end
    
    row_num = size(y1, 1);
    
    fig = figure();
    if ~isnan(options.Position)
        fig.Position = options.Position;
    end
    
    AX1 = axes(fig);
    L = gobjects(1, row_num);
    for row = 1 : 1 : row_num
        L(row) = plot(x, y1(row, :), 'LineWidth', options.LineWidth, 'Color', options.Color(row), 'Parent', AX1);
        hold on;
        plot(x, y2(row, :), '--', 'LineWidth', options.LineWidth, 'Color', options.Color(row), 'Parent', AX1);
        hold on;
    end
    
    set(AX1, 'Box', 'off');
    if ~isnan(options.XLim)
        xlim(options.XLim);
    end
    if ~isnan(options.YLim)
        ylim(options.YLim);
    end
    LGND1 = legend(AX1, L, lgnd.ColorLegend, 'location', lgnd.ColorLocation);
    hold on;
    
    AX2 = copyobj(AX1, gcf);
    delete( get(AX2, 'Children') );
    H1 = plot(NaN, 'LineWidth', options.LineWidth, 'Color', [0 0 0], 'Parent', AX2);
    H2 = plot(NaN, '--', 'LineWidth', options.LineWidth, 'Color', [0 0 0], 'Parent', AX2);

    set(AX2, 'Color', 'none', 'XTick', [], 'Visible', 'off');
    LGND2 = legend([H1 H2], lgnd.LineLegend, 'location', lgnd.LineLocation);
    LGND2_bottom = LGND1.Position(2) - LGND2.Position(4) - 0.02;
    LGND2.Position = [LGND1.Position(1) LGND2_bottom LGND2.Position(3) LGND2.Position(4)];
    LGND2.Box = 'off';
    LGND2.Color = 'none';
    hold off;
    
    xlabel(options.XLabel);
    ylabel(options.YLabel);
    title(options.Title);
end

