classdef CFGAAnalysis
    properties
        cfnl
        output_dir
    end
    
    methods
        function obj = CFGAAnalysis(cfnl, output_dir)
            obj.cfnl = cfnl;
            obj.output_dir = output_dir;
        end
        
        function generate_and_save_data(obj, N_sim)
            % Generate and save data
            t_span = obj.cfnl.nominal_trj.t; % Use the nominal trajectory time vector
            T_sim = t_span(end);
            
            % Extract system parameters
            dynamics = obj.cfnl.dynamics;
            pars = obj.cfnl.pars;
            
            % Initialize storage for simulations
            all_t_sim = cell(N_sim, 1);
            all_x_sim = cell(N_sim, 1);
            all_u_sim = cell(N_sim, 1);

            for i = 1:N_sim
                % Random initial state within bounds
                % Generate random initial state with buffer
                x0 = obj.cfnl.bnds.x_min + (obj.cfnl.bnds.x_max - obj.cfnl.bnds.x_min) .* rand(obj.cfnl.nx, 1);
                
                % Random control input function within bounds
                u = @(t) obj.cfnl.bnds.u_min + (obj.cfnl.bnds.u_max - obj.cfnl.bnds.u_min) .* rand(obj.cfnl.nu, 1);
                % Define the ODE function
                ode_func = @(t, x) dynamics(t, x, u(t), [], pars);
                
                % Simulate the system
                [t_sim, x_sim] = ode45(ode_func, t_span, x0);
                
                % Store results
                all_t_sim{i} = t_sim';
                all_x_sim{i} = x_sim';
                
                % Evaluate control input at each time step
                u_sim = arrayfun(@(t) u(t), t_sim, 'UniformOutput', false);
                u_sim = cell2mat(u_sim');
                all_u_sim{i} = u_sim;
            end
            
            % Save random simulation data
            random_filename = fullfile(obj.output_dir, 'random_simulation_data.mat');
            save(random_filename, 'all_t_sim', 'all_x_sim', 'all_u_sim');
            fprintf('Saved random simulation data to: %s\n', random_filename);
            
            % Save existing simulation data
            sim_data = obj.cfnl.plot.sim_data;

            for sim = 1:sim_data.Nsim
                t_sim = sim_data.t_sim{sim};
                x_sim = sim_data.x_sim{sim};
                u_sim = sim_data.u_sim{sim};

                % Create file name
                mat_filename = fullfile(obj.output_dir, sprintf('sim_data_%d.mat', sim));

               

                % Save the data to a MAT file
                try
                    save(mat_filename, 't_sim', 'x_sim', 'u_sim');
                catch ME
                    fprintf('Failed to save simulation data to: %s\n', mat_filename);
                    rethrow(ME);
                end
            end

            fprintf('Saved all simulation data to MAT files.\n');

            % Save nominal trajectory data separately
            nominal_trj = obj.cfnl.nominal_trj;
            nominal_t = nominal_trj.t;
            nominal_x = nominal_trj.x;
            nominal_u = nominal_trj.u;
            nominal_filename = fullfile(obj.output_dir, 'nominal_trj.mat');
            save(nominal_filename, 'nominal_t', 'nominal_x', 'nominal_u');
            fprintf('Saved nominal trajectory data to: %s\n', nominal_filename);
        end
        
        function calculate_and_save_deviations(obj)
            % Load nominal trajectory data
            load(fullfile(obj.output_dir, 'nominal_trj.mat'), 'nominal_t', 'nominal_x', 'nominal_u');

            % Load random simulation data
            load(fullfile(obj.output_dir, 'random_simulation_data.mat'), 'all_t_sim', 'all_x_sim', 'all_u_sim');

            % Number of simulations
            N_sim = length(all_t_sim);

            % Initialize storage for deviations
            state_deviations = cell(N_sim, 1);
            input_deviations = cell(N_sim, 1);

            % Compute deviations for each simulation
            for i = 1:N_sim
                t_sim = all_t_sim{i};
                x_sim = all_x_sim{i};
                u_sim = all_u_sim{i};
                
                % Interpolate nominal trajectory to match the time points of random data
                nominal_x_interp = interp1(nominal_t, nominal_x', t_sim)';
                nominal_u_interp = interp1(nominal_t, nominal_u', t_sim)';
                
                % Calculate deviations
                state_deviation = x_sim - nominal_x_interp;
                input_deviation = u_sim - nominal_u_interp;
                
                % Store deviations
                state_deviations{i} = state_deviation;
                input_deviations{i} = input_deviation;
            end

            % Save deviation data
            deviation_filename = fullfile(obj.output_dir, 'deviation_data.mat');
            save(deviation_filename, 'state_deviations', 'input_deviations');
            fprintf('Saved deviation data to: %s\n', deviation_filename);
        end
        
        function rank_condition_check(obj)
            % Load deviation data
            load(fullfile(obj.output_dir, 'deviation_data.mat'), 'state_deviations', 'input_deviations');

            % Number of simulations
            N_sim = length(state_deviations);

            % Number of states and inputs
            nx = size(state_deviations{1}, 1);
            nu = size(input_deviations{1}, 1);

            % Initialize storage for ranks
            num_time_steps = length(obj.cfnl.nominal_trj.t);
            rank_results = zeros(num_time_steps, 1);

            % Check rank at each time step
            for t = 1:num_time_steps
                % Combine state deviations and input deviations across all simulations
                state_deviation_at_t = zeros(nx, N_sim);
                input_deviation_at_t = zeros(nu, N_sim);
                
                for i = 1:N_sim
                    state_deviation_at_t(:, i) = state_deviations{i}(:, t);
                    input_deviation_at_t(:, i) = input_deviations{i}(:, t);
                end
                
                % Combine state and input deviations
                combined_deviation = [state_deviation_at_t; input_deviation_at_t];
                
                % Check the rank
                rank_results(t) = rank(combined_deviation);
            end

            % Save the rank results
            save(fullfile(obj.output_dir, 'rank_results.mat'), 'rank_results');
            fprintf('Saved rank results to: %s\n', fullfile(obj.output_dir, 'rank_results.mat'));

            % Plot the rank results
            figure;
            plot(obj.cfnl.nominal_trj.t, rank_results, 'LineWidth', 2);
            grid on;
            xlabel('Time [s]');
            ylabel('Rank');
            title('Rank of Combined State and Input Deviations');
        end
        
        function plot_data(obj)
            % Load nominal trajectory data
            load(fullfile(obj.output_dir, 'nominal_trj.mat'), 'nominal_t', 'nominal_x', 'nominal_u');

            % Load random simulation data
            load(fullfile(obj.output_dir, 'random_simulation_data.mat'), 'all_t_sim', 'all_x_sim', 'all_u_sim');

            % Number of simulations
            N_sim = length(all_t_sim);

            % Plot states
            nx = size(nominal_x, 1);
            figure;
            for j = 1:nx
                subplot(nx, 1, j);
                hold on;
                grid on;
                box on;
                
                % Plot nominal trajectory states
                plot(nominal_t, nominal_x(j, :), 'LineWidth', 2);
                
                % Plot random simulation states
                for i = 1:N_sim
                    t_sim = all_t_sim{i};
                    x_sim = all_x_sim{i};
                    plot(t_sim, x_sim(j, :), '--');
                end
                
                % Labels and titles
                title(['State ' num2str(j)]);
                xlabel('Time [s]');
                ylabel(['State ' num2str(j)]);
            end
            sgtitle('States: Random Data and Nominal Trajectory');

            % Plot control inputs
            nu = size(nominal_u, 1);
            figure;
            for k = 1:nu
                subplot(nu, 1, k);
                hold on;
                grid on;
                box on;
                
                % Plot nominal trajectory inputs
                plot(nominal_t, nominal_u(k, :), 'LineWidth', 2);
                
                % Plot random simulation inputs
                for i = 1:N_sim
                    t_sim = all_t_sim{i};
                    u_sim = all_u_sim{i};
                    plot(t_sim, u_sim(k, :), '--');
                end
                
                % Labels and titles
                title(['Input ' num2str(k)]);
                xlabel('Time [s]');
                ylabel(['Input ' num2str(k)]);
            end
            sgtitle('Inputs: Random Data and Nominal Trajectory');
        end
        
        function solve_optimization(obj)
            

            % Extract Q_max and R_max
            Q_max = obj.cfnl.fnl.Q_max;
            R_max = obj.cfnl.fnl.R_max;
            
            % Time vector for interpolation
            t_span = obj.cfnl.nominal_trj.t;
            num_time_steps = length(t_span);
            
            % Original time vector
            t_original = linspace(t_span(1), t_span(end), size(Q_max, 3));
           
            % Initialize interpolated Q_max and R_max
            Q_max_interp = zeros(obj.cfnl.nx, obj.cfnl.nx, num_time_steps);
            R_max_interp = zeros(obj.cfnl.nu, obj.cfnl.nu, num_time_steps);
            
            % Initial Regularization parameter
            epsilon = 1e-8;  % Base regularization term
            
            % Interpolate Q_max and R_max to match the time steps of the nominal trajectory
            for k = 1:num_time_steps
                Q_max_interp(:, :, k) = squeeze(interp1(t_original, permute(Q_max, [3, 1, 2]), t_span(k), 'linear', 'extrap')) + epsilon * eye(obj.cfnl.nx);
                R_max_interp(:, :, k) = squeeze(interp1(t_original, permute(R_max, [3, 1, 2]), t_span(k), 'linear', 'extrap'));
            end

            
            
            
            % Load deviation data
            load(fullfile(obj.output_dir, 'deviation_data.mat'), 'state_deviations', 'input_deviations');
            % Number of simulations
            N_sim = length(state_deviations);
            
            % Number of states and inputs
            nx = size(state_deviations{1}, 1);
            nu = size(input_deviations{1}, 1);
            
            % Number of time steps
            num_time_steps = length(obj.cfnl.nominal_trj.t);
            
            % Plot state deviations
            figure;
            for j = 1:nx
                subplot(nx, 1, j);
                hold on;
                grid on;
                box on;
                
                % Plot deviations for each simulation
                for i = 1:N_sim
                    state_deviation = state_deviations{i};
                    plot(obj.cfnl.nominal_trj.t, state_deviation(j, :), '--');
                end
                
                % Labels and titles
                title(['State Deviation ' num2str(j)]);
                xlabel('Time [s]');
                ylabel(['State Deviation ' num2str(j)]);
            end
            sgtitle('State Deviations');
            
            % Plot input deviations
            figure;
            for k = 1:nu
                subplot(nu, 1, k);
                hold on;
                grid on;
                box on;
                
                % Plot deviations for each simulation
                for i = 1:N_sim
                    input_deviation = input_deviations{i};
                    plot(obj.cfnl.nominal_trj.t, input_deviation(k, :), '--');
                end
                
                % Labels and titles
                title(['Input Deviation ' num2str(k)]);
                xlabel('Time [s]');
                ylabel(['Input Deviation ' num2str(k)]);
            end
            sgtitle('Input Deviations');

            yalmip('clear')
            % Number of simulations
            N_sim = length(state_deviations);
        
            % Number of states and inputs
            nx = size(state_deviations{1}, 1);
            nu = size(input_deviations{1}, 1);
        
            % Number of time steps
            N = length(obj.cfnl.nominal_trj.t) - 1;
        
            % Define optimization variables
            P = sdpvar(nx, nx, N+1, 'symmetric');
            Y = sdpvar(N_sim, nx, N, 'full');
            S1 = sdpvar(nx, nx, N+1, 'symmetric');
            S2 = sdpvar(nx, nx, N+1, 'symmetric');
            Constraints = [];
            


            % Adding SDP constraints for each time step
            for k = 1:N
                % Extract state and input deviations at time k
                state_deviation_at_k = zeros(nx, N_sim);
                input_deviation_at_k = zeros(nu, N_sim);
                state_deviation_at_k1 = zeros(nx, N_sim);
                for i = 1:N_sim
                    state_deviation_at_k(:, i) = state_deviations{i}(:, k);
                    input_deviation_at_k(:, i) = input_deviations{i}(:, k);
                    state_deviation_at_k1(:, i) = state_deviations{i}(:, k+1); % State deviation at k+1
                end
                disp('State deviations at each time step:');
                for i = 1:N_sim
                    disp(state_deviations{i});
                end
                
                disp('Input deviations at each time step:');
                for i = 1:N_sim
                    disp(input_deviations{i});
                end

                X_k = state_deviation_at_k;
                U_k = input_deviation_at_k;
                X_k_1 = state_deviation_at_k1; % State deviation at k+1
                 % Display the shape of X_k, U_k, and X_k_1 for debugging

                Y_k = Y(:,:,k);
                P_k = P(:,:,k);
                P_k_1 = P(:,:,k+1);
                S1_k = S1(:,:,k);
                S2_k = S2(:,:,k);
        
                % Ensure symmetry in Q_max_interp and P matrices
                Q_max_k = diag(diag(Q_max_interp(:,:,k)));
                R_max_k = diag(diag(R_max_interp(:,:,k)));
                % R_max_k = 0.1*R_max_k;
                Q_max_k2 = Q_max_k^2;
                inv_Q_max_k = inv(Q_max_k);
                Constraints = [Constraints, (Q_max_k(:) <= P_k(:))];
                Constraints = [Constraints, (P_k(:) <= 100)];
                Constraints = [Constraints, (Q_max_k(:) <= S1_k(:))];
                Constraints = [Constraints, (100 >= S1_k(:))];
                Constraints = [Constraints, (Q_max_k2(:) <= S2_k(:))];
                Constraints = [Constraints, (6000 >= S2_k(:))];
            
                
                % Original constraints
                
                schur_matrix = [P_k_1 , X_k_1 * Y_k; Y_k' * X_k_1', P_k];
                Constraints = [Constraints, schur_matrix >= 0];
                Constraints = [Constraints, X_k * Y_k == P_k];

                % Slack constraint for P_k^3 approximation
                Constraints = [Constraints, [S1_k, P_k; P_k, eye(nx)] >= 0];
                Constraints = [Constraints, [S2_k, S1_k; S1_k, P_k] >= 0];
        
                % Main inequality using S2_k instead of P^3
                funnel = [R_max_k, U_k * Y_k; Y_k' * U_k', S2_k];
                Constraints = [Constraints, funnel >= 0];

                % Constraints = [Constraints, (Q_max_k <= P_k): 'Q_max_k <= P_k'];
                % Constraints = [Constraints, (P_k <= inv(Q_max_k)): 'P_k <= inv(Q_max_k)'];


            end
            % Final constraint on P(N)
            Q_N = Q_max_interp(:,:,end);
            Q_N = diag(diag(Q_N));
            P_N = P(:, :, end);
            
            Constraints = [Constraints, (Q_N(:) <= P_N(:))];
            Constraints = [Constraints, (P(:, :, end) <= 50*eye(nx))];
        

            % Objective function using MAXDET
            objective = 0;
            for k = 1:N
                objective = objective + trace(S2(:,:,k)) + trace(S1(:,:,k));
            end
            objective = objective + 1000*P_N;

            % Solve Optimization
            options = sdpsettings('solver', 'mosek', 'verbose', 1, 'debug', 1, 'warning', 1, 'savesolveroutput', 1, 'savesolverinput', 1);
            sol = optimize(Constraints, objective, options);
            
            % Print diagnostic information
            disp('Solver output:');
            disp(sol.solveroutput);
            disp('Solver input:');
            disp(sol.solverinput);
            
            % Check the result

            if sol.problem == 0
                disp('Optimization solved successfully');
        
                % Save optimization results
                P_opt = value(P);
                Y_opt = value(Y);
                S1_opt = value(S1);
                S2_opt = value(S2);
                save(fullfile(obj.output_dir, 'optimization_solutions.mat'), 'P_opt', 'Y_opt', 'S1_opt', 'S2_opt');
        
                % Calculate K_star
                K_star = cell(1, N);
                for k = 1:N
                    % Extract the corresponding slices of Y and P
                    Y_k = Y_opt(:,:,k);
                    P_k = P_opt(:,:,k);
        
                    % Use the input deviation data for U
                    U_k = zeros(nu, N_sim);
                    for i = 1:N_sim
                        U_k(:, i) = input_deviations{i}(:, k);
                    end
        
                    % Ensure P_k is invertible
                    if rank(P_k) == size(P_k, 1)
                        % Calculate K_k
                        K_star{k} = U_k * Y_k * inv(P_k);
                    else
                        warning(['P_k is singular or nearly singular at step ', num2str(k)]);
                        % Handle the non-invertible case appropriately
                        K_star{k} = NaN; % Placeholder for singular cases
                    end
                end
        
                % Save K_star to a file
                K_star_filename = fullfile(obj.output_dir, 'K_star.mat');
                save(K_star_filename, 'K_star');
                fprintf('Saved K_star to: %s\n', K_star_filename);
        
            else
                disp(['Problem not solved, status: ' yalmiperror(sol.problem)]);
                diagnostics = sol.info;
                disp(diagnostics);
            end
end

        function plot_optimization_solutions(obj)
            % Load optimization solutions
            load(fullfile(obj.output_dir, 'optimization_solutions.mat'), 'P_opt', 'Y_opt', 'S1_opt', 'S2_opt');
            load(fullfile(obj.output_dir, 'K_star.mat'), 'K_star');
        
            % Time vector for plotting
            t_span = obj.cfnl.nominal_trj.t;
            num_time_steps = length(t_span);
        
            % Number of time steps
            N = num_time_steps - 1;
        
            % Initialize storage for norms
            P_norms = zeros(1, num_time_steps);
            S1_norms = zeros(1, num_time_steps);
            S2_norms = zeros(1, num_time_steps);
            Y_norms = zeros(1, num_time_steps - 1);
            K_norms = zeros(1, num_time_steps - 1);
        
            % Calculate norms
            for k = 1:num_time_steps
                P_norms(k) = norm(P_opt(:,:,k), 'fro');
                S1_norms(k) = norm(S1_opt(:,:,k), 'fro');
                S2_norms(k) = norm(S2_opt(:,:,k), 'fro');
                if k <= num_time_steps - 1
                    Y_norms(k) = norm(Y_opt(:,:,k), 'fro');
                    if ~isnan(K_star{k})
                        K_norms(k) = norm(K_star{k}, 'fro');
                    else
                        K_norms(k) = NaN;
                    end
                end
            end
        
            % Plot norms through time
            figure;
            subplot(5, 1, 1);
            plot(t_span, P_norms, 'LineWidth', 2);
            grid on;
            xlabel('Time [s]');
            ylabel('|P|_F');
            title('Norm of P through time');
        
            subplot(5, 1, 2);
            plot(t_span, S1_norms, 'LineWidth', 2);
            grid on;
            xlabel('Time [s]');
            ylabel('|S1|_F');
            title('Norm of S1 through time');
        
            subplot(5, 1, 3);
            plot(t_span, S2_norms, 'LineWidth', 2);
            grid on;
            xlabel('Time [s]');
            ylabel('|S2|_F');
            title('Norm of S2 through time');
        
            subplot(5, 1, 4);
            plot(t_span(1:end-1), Y_norms, 'LineWidth', 2);
            grid on;
            xlabel('Time [s]');
            ylabel('|Y|_F');
            title('Norm of Y through time');
        
            subplot(5, 1, 5);
            plot(t_span(1:end-1), K_norms, 'LineWidth', 2);
            grid on;
            xlabel('Time [s]');
            ylabel('|K^*|_F');
            title('Norm of K* through time');
        
            sgtitle('Optimization Solution Norms through Time');
        end

        function simulate_closed_loop_with_feedback(obj)
            % Load K_star
            load(fullfile(obj.output_dir, 'K_star.mat'), 'K_star');
            
            % Load nominal trajectory data
            load(fullfile(obj.output_dir, 'nominal_trj.mat'), 'nominal_t', 'nominal_x', 'nominal_u');
            
            % Extract system parameters
            dynamics = obj.cfnl.dynamics;
            pars = obj.cfnl.pars;
            
            % Time span for simulation
            t_span = nominal_t;
            num_time_steps = length(t_span);
            
            % Initial state
            x0 = [
                2.00000000000000;
                15;
                0.999999999999995;
                4.44089209850063e-15;
                1.75415237890775e-14;
                -3.77475828372553e-15
            ];
            
            % Initialize state
            x = x0;
            
            % Initialize arrays to store state and control inputs
            x_sim = zeros(obj.cfnl.nx, num_time_steps);
            u_sim = zeros(obj.cfnl.nu, num_time_steps);
            x_sim(:, 1) = x;
            
            % Simulate the system
            for k = 1:num_time_steps - 1
                
                
                % Calculate control input using K_star
                if ~isnan(K_star{k})
                    u = K_star{k} * x;
                else
                    u = nominal_u(:, k); % Use nominal input if K_star is NaN
                end
                
                % Store control input
                u_sim(:, k) = u;
                
                % Define the ODE function
                ode_func = @(t, x) dynamics(t, x, u, [], pars);
                
                % Simulate the next state
                [~, x_next] = ode45(ode_func, [t_span(k) t_span(k + 1)], x);
                
                % Update state
                x = x_next(end, :)';    
                x_sim(:, k + 1) = x;
            end
            
            % Store results
            t_sim = t_span';
            
            % Save closed-loop simulation data
            closed_loop_filename = fullfile(obj.output_dir, 'closed_loop_simulation_data.mat');
            save(closed_loop_filename, 't_sim', 'x_sim', 'u_sim');
            fprintf('Saved closed-loop simulation data to: %s\n', closed_loop_filename);
        
            % Plotting the state trajectories
            nx = size(nominal_x, 1);
            figure;
            for j = 1:nx
                subplot(nx, 1, j);
                hold on;
                grid on;
                box on;
                
                % Plot nominal trajectory states
                plot(nominal_t, nominal_x(j, :), 'LineWidth', 2, 'DisplayName', 'Nominal Trajectory');
                
                % Plot closed-loop simulation states
                plot(t_sim, x_sim(j, :), '--', 'DisplayName', 'Closed-Loop Simulation');
                
                % Labels and titles
                title(['State ' num2str(j)]);
                xlabel('Time [s]');
                ylabel(['State ' num2str(j)]);
                legend('show');
            end
            sgtitle('States: Closed-Loop Simulation and Nominal Trajectory');
        
            % Plotting the control inputs
            nu = size(nominal_u, 1);
            figure;
            for k = 1:nu
                subplot(nu, 1, k);
                hold on;
                grid on;
                box on;
                
                % Plot nominal trajectory inputs
                plot(nominal_t, nominal_u(k, :), 'LineWidth', 2, 'DisplayName', 'Nominal Trajectory');
                
                % Plot closed-loop simulation inputs
                plot(t_sim, u_sim(k, :), '--', 'DisplayName', 'Closed-Loop Simulation');
                
                % Labels and titles
                title(['Input ' num2str(k)]);
                xlabel('Time [s]');
                ylabel(['Input ' num2str(k)]);
                legend('show');
            end
            sgtitle('Inputs: Closed-Loop Simulation and Nominal Trajectory');
        end

function plot_funnels_around_nominal(obj)
    % Load K_star and optimization solutions
    load(fullfile(obj.output_dir, 'K_star.mat'), 'K_star');
    load(fullfile(obj.output_dir, 'optimization_solutions.mat'), 'P_opt');
    
    % Load nominal trajectory data
    load(fullfile(obj.output_dir, 'nominal_trj.mat'), 'nominal_t', 'nominal_x', 'nominal_u');
    
    % Define the indices for the states of interest
    states_of_interest = [1, 2, 5]; % Indices for r_{\mathcal{I},x}, r_{\mathcal{I},y}, and \theta
    state_labels = {'r_{\mathcal{I},x} [m]', 'r_{\mathcal{I},y} [m]', '\theta [deg]'};
    
    % Number of points for the ellipsoid
    num_points = 100;
    theta = linspace(0, 2*pi, num_points);
    
    % Initialize plot handles
    hNominal = [];
    hFunnel = [];
    
    % Colors
    colorNomTraj = [0, 0, 1];
    colorFunnel = [1, 0, 0];

    % Plot the funnels around the nominal trajectories for the states of interest
    figure;
    for state_idx = 1:length(states_of_interest)
        subplot(length(states_of_interest), 1, state_idx);
        hold on;

        % Plot the nominal trajectory for the current state
        hNominal = plot(nominal_t, nominal_x(states_of_interest(state_idx), :), '-', 'LineWidth', 2, 'Color', colorNomTraj);

        % Overlay the funnel as shaded area
        funnel_lower = zeros(1, length(nominal_t));
        funnel_upper = zeros(1, length(nominal_t));
        for k = 1:length(nominal_t)
            P_k = P_opt(:, :, k);
            Q_k = inv(P_k);
            funnel_radius = sqrt(Q_k(states_of_interest(state_idx), states_of_interest(state_idx)));
            
            % Calculate upper and lower bounds for the funnel
            funnel_upper(k) = nominal_x(states_of_interest(state_idx), k) + funnel_radius;
            funnel_lower(k) = nominal_x(states_of_interest(state_idx), k) - funnel_radius;
        end

        % Plot the funnel as a shaded area
        hFunnel = fill([nominal_t, fliplr(nominal_t)], [funnel_upper, fliplr(funnel_lower)], colorFunnel, 'FaceAlpha', 0.4, 'EdgeColor', 'none');

        xlabel('Time [s]');
        ylabel(state_labels{state_idx});
        title(['State ', state_labels{state_idx}, ' Trajectories']);
        grid on;
        hold off;
    end

    sgtitle('State Trajectories with Funnel');
    legend([hNominal(1), hFunnel(1)], {'Nominal', 'Funnel'}, 'Location', 'best');

    % Plot nominal control inputs
    figure;
    plot(nominal_t, nominal_u, 'k-', 'LineWidth', 2); % Nominal control inputs in black
    xlabel('Time [s]');
    ylabel('Control input values');
    title('Nominal Control Inputs');
    grid on;
end




    end
end
