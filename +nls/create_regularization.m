function func = create_regularization(measured_values, dofmap, ...
                                      m_ref, solver_type)
    % Prepare function for solving H1-regularized Gauss-Newton systems
    %
    % INPUT PARAMETER
    %   measured_values ... vector [num_observations, 1]
    %   dofmap          ... struct, representing L^2-conforming parameter space
    %   m_ref           ... vector [dofmap.dim, 1], reference parameter
    %   solver_type     ... character vector; options:
    %                           'direct'
    %                           'krylov-full-woodbury'
    %                           'krylov' (alias for 'krylov-full-woodbury')
    %                           'krylov-no-woodbury'
    %
    % OUTPUT PARAMETER
    %   func ... function of signature:
    %                [dm, iter, t1, t2, t3] = solve(d, J, m, beta);
    %
    %       dm     ... vector [dofmap.dim, 1], parameter update
    %       iter   ... number of solver iterations performed
    %       t1     ... time spent in forming Woodbury correction (only in krylov)
    %       t2     ... time spent in forming capacitance matrix (only in krylov)
    %       t3     ... time spent in Cholesky factorization of correction (only in krylov)
    %       d      ... vector [num_observations, 1], current model response
    %       J      ... matrix [num_observations, dofmap.dim], current Jacobian
    %       m      ... vector [dofmap.dim, 1], current parameter value
    %       beta   ... positive scalar, regularization parameter

    [M, D] = assemble_regularization(dofmap);
    switch solver_type
    case 'direct'
        func = create_solver_direct(M, D, m_ref, measured_values);
        return
    case {'krylov', 'krylov-full-woodbury'}
        assemble_rhs = create_rhs_assembler_krylov(measured_values, D, m_ref);
        solver = create_solver_krylov_woodbury(M, D, 'full-woodbury');
    case 'krylov-no-woodbury'
        assemble_rhs = create_rhs_assembler_krylov(measured_values, D, m_ref);
        solver = create_solver_krylov_woodbury(M, D, 'no-woodbury');
    otherwise
        error('unknown solver type "%s"', solver_type);
    end

    function [dm, iter, t1, t2, t3] = solve(d, J, m, beta)
        rhs = assemble_rhs(J, d, m, beta);
        [dm, iter, t1, t2, t3] = solver(J, beta, rhs);
    end

    func = @solve;
end


function [M, D] = assemble_regularization(dofmap_l2)
    fprintf('Assembling regularization: ');
    t = tic();

    % Build Hdiv dofmap
    mesh = dofmap_l2.mesh;
    order = dofmap_l2.element.order + 1;
    element_hdiv = fe.create_raviart_thomas_element(mesh.dim, order);

    w = warning('off', 'Mesh:ExtraConnStored');  % Mute warning
    for d = element_hdiv.get_dof_entity_dims()
        mesh.compute_connectivity(mesh.dim, d);
        mesh.clear_connectivity(d, 0);
    end
    mesh.compute_connectivity(mesh.dim, mesh.dim-1);
    mesh.clear_connectivity(mesh.dim-1, 0);
    warning(w);  % Restore warning

    dofmap_hdiv = assembling.build_dofmap(mesh, element_hdiv);

    % Assemble Hdiv mass matrix and divergence matrix
    [M, D] = assembling.assemble_hdiv_operators(dofmap_hdiv, dofmap_l2, @(x, c) 1, 0);

    % Get zero flux dofs
    bc_dofs = build_dirichlet_dofs_hdiv0(dofmap_hdiv);

    % Apply zero flux boundary conditions
    M = assembling.apply_dirichlet_bc(M, [], bc_dofs, 0);
    D(:, bc_dofs) = 0;

    fprintf('%f seconds\n', toc(t));
end


function bc_dofs = build_dirichlet_dofs_hdiv0(dofmap)

    % Apply zero flux BC nowhere
    % FIXME: Hard-coded BC here
    boundary = @(x) false;
    value = @(x) zeros(1, dofmap.mesh.dim);

    w = warning('off', 'Mesh:ExtraConnStored');  % Mute warning
    dofmap.mesh.compute_connectivity(dofmap.mesh.dim, dofmap.mesh.dim-1);
    warning(w);  % Restore warning
    dofmap.mesh.clear_connectivity(dofmap.mesh.dim-1, 0);
    dofmap.mesh.compute_boundary_facets();

    [bc_dofs, ~] = assembling.build_dirichlet_dofs(dofmap, {boundary}, {value});
end


function func = create_rhs_assembler_krylov(measured_values, D, m_ref)

    assert(~issparse(measured_values) && iscolumn(measured_values));
    assert(issparse(D) && isreal(D));
    assert(~issparse(m_ref) && iscolumn(m_ref) && isreal(m_ref));

    [N2, N1] = size(D);
    measured_values = util.complex2real(measured_values);

    function rhs = assemble_rhs(J, d, m, beta)
        fprintf('Assembling Gauss-Newton RHS: ');
        t = tic();

        assert(~issparse(J) && ismatrix(J));
        assert(~issparse(d) && iscolumn(d));
        assert(~issparse(m) && iscolumn(m) && isreal(m));
        assert(isscalar(beta) && isreal(beta) && beta > 0);

        J = util.complex2real(J);
        d = util.complex2real(d);

        rhs = zeros(N1+N2, 1);
        rhs(1:N1)       = D.'*(m_ref-m);
        rhs(N1+1:N1+N2) = (1/beta)*(J.'*(measured_values-d));

        assert(~issparse(rhs) && isreal(rhs));

        fprintf('%f seconds\n', toc(t));
    end

    func = @assemble_rhs;
end


function func = create_solver_direct(M, D, m_ref, measured_values)
    t = tic();

    assert(issparse(M) && isreal(M));
    assert(issparse(D) && isreal(D));
    assert(~issparse(m_ref) && iscolumn(m_ref) && isreal(m_ref));
    assert(~issparse(measured_values) && iscolumn(measured_values));

    N1 = size(M, 1);
    N2 = size(D, 1);
    assert(size(M, 2) == N1);
    assert(size(D, 2) == N1);

    measured_values = util.complex2real(measured_values);

    % Assemble mixed regularization matrix
    Z = sparse(N2, N2);
    A = [
        -M, D.';
         D, Z  ;
    ];

    % Compute LDLT factors
    A = decomposition(A, 'ldl');

    function [dm, iter, t1, t2, t3] = solver(d, J, m, beta)
        t = tic();

        assert(~issparse(d) && iscolumn(d));
        assert(~issparse(J) && ismatrix(J));
        assert(~issparse(m) && iscolumn(m) && isreal(m));
        assert(isscalar(beta) && isreal(beta) && beta > 0);

        J = util.complex2real(J);
        d = util.complex2real(d);

        t_ = tic();
        H = [zeros(N1, size(J, 1)); J.'];
        H = A\H;
        H = (1/beta) * H(N1+1:end, :);
        t1 = toc(t_);
        assert(~issparse(H));

        t_ = tic();
        C = J*H;
        n = size(C, 1);
        C(1:n+1:end) = C(1:n+1:end) + 1;
        t2 = toc(t_);
        assert(~issparse(C));

        t_ = tic();
        y = C \ (J*(m - m_ref) - (d - measured_values));
        t3 = toc(t_);
        assert(~issparse(y) && iscolumn(y));

        dm = m_ref - m + H*y;
        assert(~issparse(dm) && isreal(dm));
        iter = 1;

        fprintf('Solved by Cholesky-Woodbury: %f seconds\n', toc(t));
    end

    func = @solver;

    fprintf('Prepared Cholesky-Woodbury solver: %f\n', toc(t));
end


function func = create_solver_krylov_woodbury(M, D, variant)
    t = tic();

    assert(issparse(M) && isreal(M));
    assert(issparse(D) && isreal(D));

    if nargin < 3
        variant = 'full-woodbury';
    end

    [N2, N1] = size(D);

    % Inverse diagonal of M
    Mdiag_inv = 1./diag(M);
    assert(~issparse(Mdiag_inv));

    % Preconditioner for Schur: AMG with D*inv(diag(M))*DT
    Mdiag_inv_mat = spdiags(Mdiag_inv, 0, N1, N1);
    assert(issparse(Mdiag_inv_mat));
    Sdiag = D*Mdiag_inv_mat*D.';
    assert(issparse(Sdiag));
    clear Mdiag_inv_mat;
    ctrl = hsl_mi20_control;
    ctrl.st_parameter = 0.67;
    ctrl.pre_smoothing = 1;
    ctrl.post_smoothing = 1;
    ctrl.smoother = 1;
    amg_S = solving.HSLMI20(Sdiag, ctrl);
    fprintf('HSLMI20:\n'); disp(amg_S.inform);

    % Preconditioner for Hdiv block: scaling by inverse diagonal of mass
    % TODO: Try Chebyshev+Jacobi with full M?
    function x = precondition_mass(b)
        x(:, 1) = Mdiag_inv.*b;
    end

    function [dm, iter, t1, t2, t3] = solver(J, beta, rhs)

        t = tic();

        assert(~issparse(J) && ismatrix(J));
        assert(isscalar(beta) && isreal(beta) && beta > 0);
        assert(~issparse(rhs) && iscolumn(rhs) && isreal(rhs));

        J = util.complex2real(J);

        % Form correction to Schur preconditioner using Woodbury formula
        switch variant
        case 'full-woodbury'
            [S_corr, t1, t2, t3] = create_woodbury_schur_correction_(J, amg_S, beta);
        case 'no-woodbury'
            [S_corr, t1, t2, t3] = deal(@(~) 0, 0, 0, 0);
        otherwise
            error('unknown variant "%s"', variant);
        end

        % Preconditioner for perturbed Schur
        tmp_x2 = zeros(N2, 1);
        function x2 = precondition_schur(b2)
            tmp_x2(:, 1) = amg_S.precondition(b2) + S_corr(b2);
            x2 = tmp_x2;
            assert(isreal(x2));
        end

        % Block-diagonal preconditioner for the system
        tmp = zeros(N1+N2, 1);
        function x = precondition_system(b)
            tmp(N1+1:N1+N2) = precondition_schur(b(N1+1:N1+N2));
            tmp(1:N1) = precondition_mass(b(1:N1));
            x = tmp;
            assert(isreal(x));
        end

        % System operator
        tmp_y = zeros(N1+N2, 1);
        function y = Afun(z)
            tmp_y(1:N1, 1) = D.'*z(N1+1:N1+N2) - M*z(1:N1);
            % NB: Note the parentheses J.'*(J*...) !!!
            tmp_y(N1+1:N1+N2, 1) = D*z(1:N1) + J.'*(J*z(N1+1:N1+N2))/beta;
            y = tmp_y;
            assert(isreal(y));
        end

        % Run preconditioned MINRES
        tol = 1e-7;
        maxit = 2*(N1+N2);
        [x(:), flag, relres, iter] = minres(@Afun, rhs, tol, maxit, @precondition_system);

        % Handle errors
        msg = 'MINRES returned flag %d (%s) in iteration %d of %d: relres=%g (tol=%g)';
        explain = {
            'no convergence';
            'preconditioner ill-conditioned';
            'MINRES stagnated';
            'overflow/underflow';
            'preconditioner not SPD';
        };
        if flag == 1
            warning(msg, flag, explain{flag}, iter, maxit, relres, tol);
        elseif flag > 1
            error(msg, flag, explain{flag}, iter, maxit, relres, tol);
        else
            fprintf([msg, '\n'], flag, 'converged', iter, maxit, relres, tol);
        end

        % Extract the second block
        dm(:) = x(N1+1:N1+N2);

        assert(~issparse(dm) && isreal(dm));

        fprintf('Solved by Krylov-Woodbury: %f seconds\n', toc(t));
    end

    func = @solver;

    fprintf('Prepared Krylov-Woodbury solver: %f seconds\n', toc(t));
end


function [S_corr, t1, t2, t3] = create_woodbury_schur_correction_(J, amg_S, beta)

    t_ = tic();
    H = amg_S.solve(J, 'transpose_rhs');
    H = (1/sqrt(beta)) * H;
    t1 = toc(t_);
    assert(~issparse(H));

    t_ = tic();
    % FIXME: 3 nested loops have complexity O(M*N*M)
    %        Blocking can improve this, but probably not to O(M*N)
    %        Can we do something about it?
    C = J*H;
    C = (1/sqrt(beta)) * C;
    n = size(C, 1);
    C(1:n+1:end) = C(1:n+1:end) + 1;
    t2 = toc(t_);
    assert(~issparse(C));

    t_ = tic();
    R = chol(C);
    t3 = toc(t_);
    assert(~issparse(R));

    % NB: Operator S_corr cannot be assembled! It is dense!
    %     Hence using matrix free action on vector b2.
    S_corr = @(b2) -(H*(R\(R.'\(H.'*b2))));

end
