clear; clear cfga;
% Load the .mat file from the 'data' folder
data = load('data/planarquad_scvx.mat');
% Define new values
newN = 30; % example new value for scvx.ctrl.N
newP = 6.5; % example new value for scvx.output.p, adjust dimensions as needed
% Set the new values
data.scvx.ctrl.N = newN;
data.scvx.output.p = newP;

% Save the modified data back to the .mat file in the 'data' folder
save('data/planarquad_scvx.mat', '-struct', 'data');
%%
run('utils/set_path.m')

% instantiate a CFGA object
cfnl = cfga('planarquad');

% general options
cfnl.opts.cvrg_min      = 3;
cfnl.opts.max_iter      = 100;
cfnl.opts.decay_rate    = 0.01;
cfnl.opts.lmi_tol       = 0;
cfnl.plot.make_rss      = true;

% variable contraction parameters
cfnl.opts.contract_min      = 0.7;
cfnl.opts.contract_width    = 15;

% convergence tolerances
cfnl.opts.cvrg_tol_a = deg2rad(3);
cfnl.opts.cvrg_tol_w = deg2rad(3);

% attach dynamics and linearization functions
cfnl.dynamics  = @dynamics;
cfnl.linearize = @linearize;

% set SCvx data from stored trajectory
cfnl.set_scvx_data();

% synthesize funnels
cfnl.synthesize_funnel(20);

% simulate some test cases ( # cases, t_{start} )
cfnl.get_sim_data(25,0);

% make the desired plots
cfnl.plot.make_plots(cfnl);
%%
clear CFGAAnalysis;
% Create CFGAAnalysis object
analysis = CFGAAnalysis(cfnl, 'simulation_data');

% Generate and save data
analysis.generate_and_save_data(10);
%%
% Calculate and save deviations
analysis.calculate_and_save_deviations();

% Perform rank condition check
analysis.rank_condition_check();
% Plot the data
analysis.plot_data();
% Solve the optimization problem
analysis.solve_optimization();
%%
analysis.plot_optimization_solutions();
analysis.simulate_closed_loop_with_feedback();
analysis.plot_funnels_around_nominal();