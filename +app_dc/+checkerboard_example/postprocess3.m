%compose_resistivity_figure_2d('dir', [4, 8, 16, 32, 64, 128, 256]);
%compose_resistivity_figure_2d('fw',  [4, 8, 16, 32, 64, 128, 256]);
%compose_resistivity_figure_2d('nw',  [4, 8, 16, 32, 64]);
compose_resistivity_figure_2d('dir', [4, 8, 16, 32]);
compose_resistivity_figure_2d('fw',  [4, 8, 16, 32, 64]);
compose_resistivity_figure_2d('nw',  [4, 8]);


function compose_resistivity_figure_2d(tag, ns)

    tag = sprintf('checkerboard-resistivity-2d%s', tag);
    infile = sprintf('%s-n%%04d.fig', tag);
    outfile = sprintf('%s.pdf', tag);

    aspect_ratio = numel(ns)/2 * 0.5;
    hspace = 0.02;
    cbar_height = 0.08;
    vspace = hspace / aspect_ratio;
    height = (1-cbar_height)/numel(ns);
    dpi = 92;
    paperw = 5.125;
    marker_sizes = dpi/92 * [12, 8, 3, 1, 0.5, 0, 0, 0, 0];

    fig = figure();
    axs = gobjects(0);
    [lo, hi] = deal(+inf, -inf);

    cm_true = parula();
    cm_true(end, :) = [1, 0, 0];  % overshoot is red
    cm_inv = parula();

    for row = 1:numel(ns)
        f = openfig(sprintf(infile, ns(row)), 'invisible');
        figure(fig);

        pos = [0, 1-cbar_height-row*height+vspace/2, 0.5-hspace, height-vspace];
        axes('Units', 'normalized', 'Position', pos);
        axis('image', 'off');
        copyobj(allchild(f.Children(6)), gca);
        colormap(gca, cm_true);
        axs(end+1) = gca;

        num_electrodes_1d = 1+4*ns(row);
        electrode_coords = create_electrode_coords_1d(num_electrodes_1d);
        marker_size = marker_sizes(row);
        if marker_size > 0
            hold on;
            plot(electrode_coords(1, :), electrode_coords(2, :), 'k|', 'MarkerSize', marker_size);
            hold off;
        end

        pos = [0.5+hspace/2, 1-cbar_height-row*height+vspace/2, 0.5-hspace, height-vspace];
        axes('Units', 'normalized', 'Position', pos);
        axis('image', 'off');
        copyobj(allchild(f.Children(2)), gca);
        colormap(gca, cm_inv);
        axs(end+1) = gca;

        limits = caxis;
        lo = min(lo, limits(1));
        hi = max(hi, limits(2));

        close(f);
    end

    for ax = axs
        caxis(ax, [lo, hi]);
    end

    cbar = colorbar('Location', 'northoutside', 'TickLabelInterpreter', 'latex');
    cbar.FontName = 'cmr';
    cbar.FontSize = 12;
    cbar.TickLabelInterpreter = 'latex';
    %cbar.Box = 'off';
    %cbar.LineWidth = 1e-323;
    %cbar.TickLength = 0;
    cbar.Position = [0.1, 1-0.7*cbar_height, 0.75, 0.2*cbar_height];
    anno_position = [0.86, 1-0.7*cbar_height, 0.2*cbar_height*aspect_ratio, 0.2*cbar_height];
    line_x = [0.86 + 0.1*cbar_height*aspect_ratio, 0.86 + 0.1*cbar_height*aspect_ratio];
    line_y = [1-0.7*cbar_height - 0.01, 1-0.7*cbar_height - 0.001];
    annotation(fig, 'rectangle', anno_position, 'FaceColor', 'r', 'LineStyle', 'none');
    annotation('textarrow', line_x, line_y, 'String', '7000', 'HeadStyle', 'none', ...
               'FontName', 'cmr', 'FontSize', 12, 'Interpreter', 'latex', 'TextMargin', 0.05);

    set(fig, 'PaperUnits', 'points');
    set(fig, 'PaperSize', dpi*[paperw, paperw*aspect_ratio]);
    set(fig, 'PaperPosition', dpi*[0, 0, paperw, paperw*aspect_ratio]);
    print(fig, outfile, '-dpdf');
    close(fig);

end


function coords = create_electrode_coords_1d(num_electrodes)
    coords = zeros(2, num_electrodes);
    coords(1, :) = create_electrode_coords_(num_electrodes);
end


function coords = create_electrode_coords_(num_electrodes)
    coords = linspace(-50, 50, num_electrodes);
end
