classdef HSLMI20 < handle

    properties (GetAccess = public, SetAccess = private)
        inform
        handle
    end

    methods (Access = public)

        function obj = HSLMI20(A, control)
            % Construct HSLMI20 AMG solver
            %
            % INPUT ARGS
            %   A ... sparse square matrix
            % OPTIONAL ARGS
            %   control ... HSLMI20 parameters obtained with hsl_mi20_control()

            if ~obj.is_installed()
                error('HSL_MI20 package not installed!');
            end

            hsl_mi20_startup();

            if nargin == 1
                control = hsl_mi20_control();
            end

            [obj.inform, obj.handle] = hsl_mi20('setup', A, control);
        end

        function delete(obj)
            fprintf('Deleting hsl_mi20 handle %d\n', obj.handle);
            hsl_mi20('destroy', obj.handle);
        end

        function x = solve(obj, b, varargin)
            % Solve equation A*x = b
            %
            % b can have multiple columns;  optional argument
            % 'transpose_rhs' can be given to transpose b

            % Parse optional arguments
            if nargin == 2
                transpose_rhs = false;
            elseif nargin == 3
                switch varargin{1}
                case 'transpose_rhs'
                    transpose_rhs = true;
                case default
                    error('Argument ''%s'' not understood', varargin{1});
                end
            else
                error('Unexpected number of input arguments');
            end

            % Allocate return value and solve column-by-column
            if transpose_rhs
                x = zeros(size(b, 2), size(b, 1));
                for j = 1:size(b, 1)
                    [x(:, j), ~] = hsl_mi20('precondition', b(j, :), obj.handle);
                end
            else
                x = zeros(size(b));
                for j = 1:size(b, 2)
                    [x(:, j), ~] = hsl_mi20('precondition', b(:, j), obj.handle);
                end
            end
        end

        function x = precondition(obj, b)
            % Solve equation A*x = b
            %
            % This method does not support multiple right-hand sides.
            % Compared to the 'solve' method, this directly executes
            % hsl_mi20('precondition', ...) without extra overhead.
            [x, ~] = hsl_mi20('precondition', b, obj.handle);
        end

    end

    methods (Static)

        function result = is_installed()
            try
                hsl_mi20_startup();
                result = true;
            catch ME
                if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
                    result = false;
                else
                    rethrow(ME);
                end
            end
        end

    end

end
