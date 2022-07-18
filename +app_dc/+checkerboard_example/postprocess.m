write_table('2ddir');
write_table('2dfw');
write_table('2dnw');
write_table('3ddir');
write_table('3dfw');
write_table('3dnw');


function write_table(tag)

    matfile = load(sprintf('checkerboard-timings-%s.mat', tag));
    timings = matfile.timings;

    t_observe = vertcat(timings.t_observe);  %#ok<NASGU>
    t_normal = vertcat(timings.t_normal);
    t_woodbury1 = vertcat(timings.t_woodbury1);
    t_woodbury2 = vertcat(timings.t_woodbury2);
    t_woodbury3 = vertcat(timings.t_woodbury3);
    iter_normal = vertcat(timings.iter_normal);
    num_dofs_fw = vertcat(timings.num_dofs_fw);  %#ok<NASGU>
    num_dofs_inv = vertcat(timings.num_dofs_inv);
    num_obs = vertcat(timings.num_obs);
    num_electrodes = vertcat(timings.num_electrodes);
    misfit = vertcat(timings.misfit);  %#ok<NASGU>
    beta = vertcat(timings.beta);  %#ok<NASGU>
    n = vertcat(timings.n);

    [num_rows, num_cols] = size(t_normal);
    if size(n, 2) == 1
        n = [n, ones(size(n))];
    end
    dup = @(arr) repelem(arr, num_cols, 1);
    iter = repmat(transpose(1:num_cols), num_rows, 1);
    flat = @(arr) reshape(transpose(arr), [], 1);
    table = [dup(n), dup(num_electrodes), dup(num_dofs_inv), dup(num_obs), ...
             iter, ...
             flat(t_woodbury1), flat(t_woodbury2), flat(t_woodbury3), ...
             flat(iter_normal), flat(t_normal)];
    filename = sprintf('checkerboard-timings-flat-%s.tex', tag);
    write_latex(filename, table);
    fprintf('Saved ''%s''\n', filename);

end


function write_latex(filename, arr)

    % Generate latex code
    prev = sympref('FloatingPointOutput', true);
    code = latex(sym(arr));
    sympref('FloatingPointOutput', prev);

    % Strip all prolog and epilog code (keep just numbers)
    code = replace(code, '\left(', '');
    code = replace(code, '\right)', '');
    code = regexprep(code, '\\begin{array}{c*}', '');
    code = replace(code, '\end{array}', '');
    code = strip(code);

    % Write to file
    [fid, msg] = fopen(filename, 'wt');
    if fid == -1
        error('Opening "%s" for writing failed with reason:\n%s', filename, msg);
    end
    fwrite(fid, code, 'char');
    fclose(fid);

end
