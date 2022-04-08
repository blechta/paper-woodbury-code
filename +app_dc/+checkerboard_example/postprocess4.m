% Select case to plot
n = [2, 4];
Nele = 289;
N = 452736;
M = 1564;
t = sprintf('3D \\normalfont ($N^\\mathrm{ele}=%d$, $N=%d$, $M=%d$)', Nele, N, M);

plot_performance_characteristics_compare_beta('t_normal', n);
ylabel('$t_{\mathrm{norm}}$ [s]', 'interpreter', 'latex');
title(t, 'interpreter', 'latex');
saveas(gcf, 'checkerboard-3d-beta-tnormal.pdf');
export_tikz('checkerboard-3d-beta-tnormal.tex');

plot_performance_characteristics_compare_beta('iter_normal', n);
ylabel('$n_{\mathrm{iter}}$', 'interpreter', 'latex');
title(t, 'interpreter', 'latex');
saveas(gcf, 'checkerboard-3d-beta-niter.pdf');
export_tikz('checkerboard-3d-beta-niter.tex');

plot_performance_characteristics_compare_beta('misfit', n);
ylabel('misfit', 'interpreter', 'latex');
title(t, 'interpreter', 'latex');
saveas(gcf, 'checkerboard-3d-beta-misfit.pdf');
export_tikz('checkerboard-3d-beta-misfit.tex');


function plot_performance_characteristics_compare_beta(yfield, n);

    figure();
    labels = {};

    tags = {'3dfw-beta50', '3dfw-beta45', '3dfw-beta40', '3dfw-beta35', '3dfw-beta30'};
    [beta, y] = extract_data(tags, yfield, n);
    p = loglog(beta, y, 'x');
    p(1).Marker = 'x';
    p(2).Marker = '+';
    for i = 1:size(y, 2)
        labels{end+1} = sprintf('$i=%d$ (%s)', i, tag_to_desc('3dfw'));  %#ok<AGROW>
    end

    hold on;

    tags = {'3dnw-beta50', '3dnw-beta45', '3dnw-beta40', '3dnw-beta35', '3dnw-beta30'};
    [beta, y] = extract_data(tags, yfield, n);
    p = loglog(beta, y, 'x');
    p(1).Marker = 'd';
    p(2).Marker = 's';
    for i = 1:size(y, 2)
        labels{end+1} = sprintf('$i=%d$ (%s)', i, tag_to_desc('3dnw'));  %#ok<AGROW>
    end

    xlabel('$\beta$', 'interpreter', 'latex');
    legend(labels, 'Location', 'northeastoutside', 'interpreter', 'latex');

    ytickangle(90);

    ax = gca;
    ax.PlotBoxAspectRatioMode = 'manual';
    ax.PlotBoxAspectRatio = [1, 1, 1];

end


function [data_x, data_y] = extract_data(tags, y_name, n)

    data_x = double.empty(0, 2);  % hard-coded number of GN iterations!
    data_y = double.empty(0, 2);  % hard-coded number of GN iterations!

    for tag = tags
        tag = tag{:};
        matfile = load(sprintf('checkerboard-timings-%s.mat', tag));
        timings = matfile.timings;
        x = vertcat(timings.('beta'));
        y = vertcat(timings.(y_name));
        [~, ind] = ismember(n, vertcat(timings.n), 'rows');
        if ind == 0
            continue
        end
        data_x(end+1, :) = x(ind, :);
        data_y(end+1, :) = y(ind, :);
    end
end


function desc = tag_to_desc(tag)
    switch tag(3:4)
    case 'di'
        desc = 'Algorithm~\ref{alg:gn-direct}';
    case 'fw'
        desc = 'Algorithm~\ref{alg:gnstep-iterative}';
    case 'nw'
        desc = 'Algorithm~\ref{alg:gnstep-iterative-nw}';
    end
end


function export_tikz(filename)
    % Try exporting TikZ using matlab2tikz package or fail gracefully

    opts = {
        'showInfo'; false;
        'width'; '\figwidth';
        'height'; '\figheight';
        'parseStrings'; false;
        'extraTikzpictureOptions'; {'trim axis left', 'trim axis right'};
        'extraAxisOptions'; {'xlabel near ticks', 'ylabel near ticks', ...
                             'yticklabel style={anchor=south west,xshift=-2.5mm}'};
    };

    try
        matlab2tikz(filename, opts{:});
    catch ME
        switch ME.identifier
        case 'MATLAB:UndefinedFunction'
            warning('matlab2tikz package not installed; skipping TikZ export');
        otherwise
            rethrow(ME)
        end
    end
end
