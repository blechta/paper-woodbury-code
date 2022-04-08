plot_performance_characteristics_compare({'2dfw', '2dnw', '2ddir'}, 0.2);
title('2D', 'interpreter', 'latex');
saveas(gcf, 'checkerboard-2d-normal.pdf');
export_tikz('checkerboard-2d-normal.tex');

plot_performance_characteristics_compare({'3dfw', '3dnw', '3ddir'}, 0.7);
title('3D \normalfont ($\beta=10^5$)', 'interpreter', 'latex');
%plot_performance_characteristics_compare({'3dfw-beta30', '3dnw-beta30'}, 0.7);
%title('3D \normalfont ($\beta=10^3$)', 'interpreter', 'latex');
saveas(gcf, 'checkerboard-3d-normal.pdf');
export_tikz('checkerboard-3d-normal.tex');


function plot_performance_characteristics_compare(tags, shift)

    figure();
    labels = {};
    markers_seq = 'x+dsv^';

    for tag = tags
        tag = tag{:};
        matfile = load(sprintf('checkerboard-timings-%s.mat', tag));
        timings = matfile.timings;
        t_normal = vertcat(timings.t_normal);
        num_dofs_inv = vertcat(timings.num_dofs_inv);
        num_obs = vertcat(timings.num_obs);

        p = loglog(num_obs.*num_dofs_inv, t_normal(:, :), 'x');
        hold on;
        p(1).Marker = markers_seq(1);
        p(2).Marker = markers_seq(2);
        markers_seq = circshift(markers_seq, -2);

        for i = 1:size(t_normal, 2)
            labels{end+1} = sprintf('$i=%d$ (%s)', i, tag_to_desc(tag));  %#ok<AGROW>
        end

    end

    tag = tags{1};
    matfile = load(sprintf('checkerboard-timings-%s.mat', tag));
    timings = matfile.timings;
    t_normal = vertcat(timings.t_normal);
    num_dofs_inv = vertcat(timings.num_dofs_inv);
    num_obs = vertcat(timings.num_obs);

    [slope_x_, slope_y_] = get_slope(num_obs.*num_dofs_inv, t_normal(:, end), 1);
    slope_y_ = shift*slope_y_;
    loglog(slope_x_, slope_y_, '--');

    x2 = num_dofs_inv.*num_obs;
    y2 = num_dofs_inv.*num_obs.^2;
    y2 = y2 * slope_y_(end) / y2(end);
    ymin = min(t_normal(:));
    x2 = x2(y2>=ymin);
    y2 = y2(y2>=ymin);
    loglog(x2, y2, '-.');

    labels{end+1} = 'slope $O(M N)$';
    labels{end+1} = 'slope $O(M^2 N)$';
    legend(labels, 'Location', 'northeastoutside', 'interpreter', 'latex');

    xlabel('$MN$', 'interpreter', 'latex');
    ylabel('$t_{\mathrm{norm}}$ [s]', 'interpreter', 'latex');

    ytickangle(90);

    ax = gca;
    ax.PlotBoxAspectRatioMode = 'manual';
    ax.PlotBoxAspectRatio = [1, 1, 1];

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


function [x_, y_] = get_slope(x, y, k)
    x_ = [min(x(:)), max(x(:))];
    y_ = x_.^k;
    y_ = y_ * max(y(:))/max(y_(:));
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
