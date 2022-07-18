ns = [4, 8, 16, 32, 64];

figure();
pos = get(gcf(), 'Position');
set(gcf(), 'Position', [pos(1), pos(2), 1000, 600]);

layout = tiledlayout(1, numel(ns));
%layout.Padding = 'compact';
%layout.TileSpacing = 'none';

for i = 1:numel(ns)

    n = ns(i);

    nexttile();
    axis('tight');

    num_electrodes_1d = 1+4*n;
    [Mtx, Mrx] = app_dc.create_electrode_configuration('pdp', [2, 4, 8], num_electrodes_1d);

    num_configs = size(Mtx, 2);
    marker_size = max(0.5, min(5, 410/num_configs));

    tstr1 = sprintf('$N^{\\mathrm{ele}} = %d$', num_electrodes_1d);
    tstr2 = sprintf('$M = %d$', num_configs);
    t = title({tstr1, tstr2}, 'rotation', 0, 'interpreter', 'latex');
    tpos = get(t, 'position');
    set(t, 'position', tpos+[0.5*num_electrodes_1d, -0.02*num_configs, 0]);

    hold('on');
    if i == 1 || i == 2
        unused = Mtx == 0 & Mrx == 0;
        spy(unused.', 'k.', 1);
    end
    spy(Mtx.', 'rv', 0.7*marker_size);
    spy(Mrx.', 'bx', marker_size);
    hold('off');

    if i == 1
        ylabel('$i$', 'interpreter', 'latex');
    end

    xlim([1, num_electrodes_1d]);
    ylim([1, inf]);

    [a, b] = deal(1, num_electrodes_1d);
    xticks([a, 0.75*a+0.25*b, 0.5*a+0.5*b, 0.25*a+0.75*b, b]);
    xticklabels({'-50', '-25', '0', '25', '50'});
    xtickangle(0);
    ytickangle(90);

    xlabel('$x$', 'rotation', 0, 'interpreter', 'latex');

    set(gca, 'TickLabelInterpreter', 'latex');
    set(gca, 'XDir', 'normal');
    set(gca, 'YDir', 'reverse');

end

set(gcf(), 'PaperOrientation', 'landscape');
set(gcf(), 'PaperUnits', 'normalized');
set(gcf(), 'PaperPosition', [0 0 1 1]);
export_pdf('checkerboard-survey.pdf');
export_tikz('checkerboard-survey.tex');


function export_pdf(filename)
    w = warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode');
    saveas(gcf, filename);
    warning(w);
    fprintf('Saved ''%s''\n', filename);
end


function export_tikz(filename)
    % Try exporting TikZ using matlab2tikz package or fail gracefully

    try
        matlab2tikz(filename, 'showInfo', false);
    catch ME
        switch ME.identifier
        case 'MATLAB:UndefinedFunction'
            warning('matlab2tikz package not installed; skipping TikZ export');
            return
        otherwise
            rethrow(ME)
        end
    end
    fprintf('Saved ''%s''\n', filename);
end
