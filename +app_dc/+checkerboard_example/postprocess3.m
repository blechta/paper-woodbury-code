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
    cbar.Position = [0.05, 1-0.7*cbar_height, 0.9, 0.2*cbar_height];

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
