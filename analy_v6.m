%% Part 1: 主计算程序
% 扫描v123参数并计算距离，保存结果
% ========================================================================

clear; clc;

% 添加路径并加载鞍点数据
addpath('./Tools/');
load('spall_ODE_siderophore_iron_partition.mat');

% 基本参数设置
R_supply = 1;
a_11 = 0.3; a_22 = 0.4; a_33 = 0.5;
alpha_matrix = [a_11, 0, 0; 0, a_22, 0; 0, 0, a_33];
v321 = 0;

N_spe = 3;
N_sid = 3;

% 参数扫描设置
v123_list = 0.64 + 0.005 * (1:32);
n_params = length(v123_list);

% 指定需要详细绘图的参数值
plot_v123_targets = [0.7150, 0.7500, 0.8];

% 时间设置
t_relax = 5000;
t_record = 2000;

% 震荡判断阈值
oscillation_var_threshold = 1e-4;

% 多初值设置
n_initial_conditions = 5;

% 鞍点筛选参数
Fe_threshold = 0.1;

% 存储结果
results = struct();
results.v123_list = v123_list;
results.param_results = cell(n_params, 1);
results.plot_v123_targets = plot_v123_targets;

% 主循环
fprintf('开始参数扫描...\n');
fprintf('总共 %d 个参数点\n', n_params);
fprintf('每个参数尝试 %d 个不同初值\n', n_initial_conditions);
fprintf('鞍点筛选条件: Fe < %.2f\n', Fe_threshold);
fprintf('震荡判断阈值: 方差 > %.2e\n', oscillation_var_threshold);
fprintf('================================================\n');

for i = 1:n_params
    fprintf('\n[%d/%d] v123 = %.4f\n', i, n_params, v123_list(i));
    
    % 筛选鞍点
    current_all_sps = spall{i}{2}; 
    n_sps = size(current_all_sps, 2);
    
    valid_sps = cell(0, 1); 
    valid_sp_indices = []; 
    
    for j = 1:n_sps
        sp = current_all_sps(:, j);
        Fe_value = sp(7);
        
        if Fe_value < Fe_threshold
            valid_sps{end+1} = sp; 
            valid_sp_indices(end+1) = j; 
        end
    end
    
    n_valid_sps = length(valid_sps);
    fprintf('  原始一阶鞍点数: %d, 筛选后: %d (Fe < %.2f)\n', ...
            n_sps, n_valid_sps, Fe_threshold);
    
    if n_valid_sps == 0
        warning('  警告: 没有满足条件的鞍点！');
        results.param_results{i}.valid = false;
        results.param_results{i}.has_oscillation = false;
        continue;
    end
    
    fprintf('  有效鞍点的Fe值: ');
    for j = 1:n_valid_sps
        fprintf('%.4f ', valid_sps{j}(7));
    end
    fprintf('\n');
    
    % 构建ODE系统
    v123 = v123_list(i);
    vself = 1 - v123 - v321;
    v_matrix = [vself, v123, v321;
                v321, vself, v123;
                v123, v321, vself];
    opt.alpha_sid_pro = alpha_matrix;
    opt.v_sid_rec = v_matrix;
    opt.e = 10;
    opt.u = 1;
    opt.migr = 0;
    opt.gamma = 1;
    opt.d = 0.1;
    opt.R_sup = R_supply;
    f = @(t, y) ODE_siderophore_iron_partition(y, opt);
    
    % 尝试多个初值
    fprintf('  尝试 %d 个不同初值...\n', n_initial_conditions);
    
    oscillating_trajectories = cell(0, 1);
    
    for init_idx = 1:n_initial_conditions
        rng(i * 100 + init_idx);
        
        biomass_init = 0.05 + rand(N_spe, 1) * 0.9;
        iron_init = 0.05 + rand(N_sid, 1) * 0.9;
        resource_init = 0.5 + rand * 0.5;
        y_init = [biomass_init; iron_init; resource_init];
        
        options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, ...
                         'NonNegative', 1:(N_spe+N_sid+1));
        try
            [~, z_relax] = ode15s(f, [0, t_relax], y_init, options);
            y0_record = z_relax(end, :)';
            
            [t_traj, z_traj] = ode15s(f, [0, t_record], y0_record, options);
            
            half_idx = floor(size(z_traj, 1) / 2);
            biomass_traj = z_traj(half_idx:end, 1:N_spe);
            biomass_var = var(biomass_traj(:));
            
            if biomass_var > oscillation_var_threshold
                traj_data = struct();
                traj_data.t = t_traj;
                traj_data.z = z_traj;
                traj_data.biomass_var = biomass_var;
                traj_data.init_idx = init_idx;
                oscillating_trajectories{end+1} = traj_data;
                fprintf('    初值 %d: 发现震荡 (方差=%.4e)\n', init_idx, biomass_var);
            else
                fprintf('    初值 %d: 无震荡 (方差=%.4e)\n', init_idx, biomass_var);
            end
        catch ME
            warning('    初值 %d: ODE求解失败 - %s', init_idx, ME.message);
        end
    end
    
    n_osc_found = length(oscillating_trajectories);
    fprintf('  --> 找到 %d 条震荡轨迹\n', n_osc_found);
    
    if n_osc_found == 0
        results.param_results{i}.valid = true;
        results.param_results{i}.has_oscillation = false;
        results.param_results{i}.v123 = v123;
        results.param_results{i}.n_valid_sps = n_valid_sps;
        continue;
    end
    
    % 选择方差最大的轨迹
    max_var = -inf;
    best_traj_idx = 1;
    for traj_idx = 1:n_osc_found
        if oscillating_trajectories{traj_idx}.biomass_var > max_var
            max_var = oscillating_trajectories{traj_idx}.biomass_var;
            best_traj_idx = traj_idx;
        end
    end
    
    selected_traj = oscillating_trajectories{best_traj_idx};
    t_traj = selected_traj.t;
    z_traj = selected_traj.z;
    biomass_var = selected_traj.biomass_var;
    
    fprintf('  选择初值 %d 的轨迹进行分析 (方差=%.4e)\n', ...
            selected_traj.init_idx, biomass_var);
    
    % 计算两种距离
    sp_results = struct();
    sp_results.n_valid_sps = n_valid_sps;
    sp_results.saddle_points = valid_sps;
    sp_results.sp_indices = valid_sp_indices;
    sp_results.biomass_var = biomass_var;
    
    sp_results.min_distances_full = zeros(n_valid_sps, 1);
    sp_results.min_distances_biomass = zeros(n_valid_sps, 1);
    sp_results.closest_points_full = cell(n_valid_sps, 1);
    sp_results.closest_points_biomass = cell(n_valid_sps, 1);
    sp_results.closest_times_full = zeros(n_valid_sps, 1);
    sp_results.closest_times_biomass = zeros(n_valid_sps, 1);
    
    fprintf('  计算到各鞍点的距离:\n');
    
    for j = 1:n_valid_sps
        sp = valid_sps{j};
        n_points = size(z_traj, 1);
        
        % 全空间距离
        distances_full = zeros(n_points, 1);
        for k = 1:n_points
            distances_full(k) = norm(z_traj(k, :)' - sp);
        end
        
        [min_dist_full, min_idx_full] = min(distances_full);
        sp_results.min_distances_full(j) = min_dist_full;
        sp_results.closest_points_full{j} = z_traj(min_idx_full, :)';
        sp_results.closest_times_full(j) = t_traj(min_idx_full);
        
        % Biomass空间距离
        distances_biomass = zeros(n_points, 1);
        for k = 1:n_points
            distances_biomass(k) = norm(z_traj(k, 1:3)' - sp(1:3));
        end
        
        [min_dist_biomass, min_idx_biomass] = min(distances_biomass);
        sp_results.min_distances_biomass(j) = min_dist_biomass;
        sp_results.closest_points_biomass{j} = z_traj(min_idx_biomass, :)';
        sp_results.closest_times_biomass(j) = t_traj(min_idx_biomass);
        
        fprintf('    鞍点 #%d: 全空间=%.6f, Biomass=%.6f\n', ...
                j, min_dist_full, min_dist_biomass);
    end
    
    [overall_min_dist_full, overall_min_sp_idx_full] = min(sp_results.min_distances_full);
    [overall_min_dist_biomass, overall_min_sp_idx_biomass] = min(sp_results.min_distances_biomass);
    
    fprintf('  --> 最短距离: 全空间=%.6f (SP#%d), Biomass=%.6f (SP#%d)\n', ...
            overall_min_dist_full, overall_min_sp_idx_full, ...
            overall_min_dist_biomass, overall_min_sp_idx_biomass);
    
    % 保存结果
    sp_results.trajectory_t = t_traj;
    sp_results.trajectory_z = z_traj;
    sp_results.v123 = v123;
    sp_results.valid = true;
    sp_results.has_oscillation = true;
    sp_results.overall_min_dist_full = overall_min_dist_full;
    sp_results.overall_min_sp_idx_full = overall_min_sp_idx_full;
    sp_results.overall_min_dist_biomass = overall_min_dist_biomass;
    sp_results.overall_min_sp_idx_biomass = overall_min_sp_idx_biomass;
    
    results.param_results{i} = sp_results;
end

fprintf('\n================================================\n');
fprintf('参数扫描完成！\n');

% 保存结果
save('saddle_distance_analysis_final.mat', 'results');
fprintf('结果已保存到: saddle_distance_analysis_final.mat\n\n');

% 输出统计信息
print_statistics(results, v123_list, n_params);

fprintf('\n提示: 可以运行 Part 2 进行绘图分析\n');

%% Part 2: 独立绘图模块
% 加载已保存的数据进行可视化，可多次运行
% ========================================================================

clear; clc;
load('saddle_distance_analysis_final.mat', 'results');

fprintf('已加载数据，开始绘图...\n');

% 提取参数
v123_list = results.v123_list;
n_params = length(v123_list);
plot_v123_targets = results.plot_v123_targets;
N_spe = 3;

% 绘制指定参数的详细图
% for i = 1:n_params
%     if isfield(results.param_results{i}, 'has_oscillation') && ...
%        results.param_results{i}.has_oscillation
%         
%         v123 = results.param_results{i}.v123;
%         
%         if any(abs(v123 - plot_v123_targets) < 1e-6)
%             fprintf('绘制 v123 = %.4f 的详细图...\n', v123);
%             sp_results = results.param_results{i};
%             plot_parameter_results(i, v123, sp_results.trajectory_z, ...
%                                  sp_results.trajectory_t, sp_results, N_spe);
%         end
%     end
% end

% 绘制三参数对比图（Biomass空间）
% plot_three_params_comparison_biomass(results, plot_v123_targets, N_spe);

% 绘制三参数对比图（全空间）
plot_three_params_comparison_full(results, plot_v123_targets, N_spe);

% 绘制汇总统计图
% plot_summary_full(results, v123_list, n_params);
% plot_summary_biomass(results, v123_list, n_params);

% 绘制柱状图对比
% plot_barplot_comparison(results, plot_v123_targets);

fprintf('绘图完成！\n');

%% 辅助函数定义
% ========================================================================

function plot_parameter_results(param_idx, v123, z_traj, t_traj, sp_results, N_spe)
    % 绘制单个参数的详细结果
    
    figure('Name', sprintf('Parameter %d v123=%.4f', param_idx, v123), ...
           'Position', [50, 50, 1400, 700]);
    
    n_valid_sps = sp_results.n_valid_sps;
    colors_sp = lines(n_valid_sps);
    
    % 子图1: 全空间距离 3D相图
    subplot(2, 3, 1);
    h_traj = plot3(z_traj(:, 1), z_traj(:, 2), z_traj(:, 3), 'b-', 'LineWidth', 1);
    hold on;
    
    h_legend = [];
    legend_labels = {};
    h_legend(end+1) = h_traj;
    legend_labels{end+1} = 'Trajectory';
    
    for j = 1:n_valid_sps
        sp = sp_results.saddle_points{j};
        cp_full = sp_results.closest_points_full{j};
        
        h_sp = plot3(sp(1), sp(2), sp(3), '*', 'Color', colors_sp(j, :), ...
                     'MarkerSize', 15, 'LineWidth', 2);
        
        plot3(cp_full(1), cp_full(2), cp_full(3), 'o', 'Color', colors_sp(j, :), ...
              'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', colors_sp(j, :));
        
        plot3([sp(1), cp_full(1)], [sp(2), cp_full(2)], [sp(3), cp_full(3)], ...
              '--', 'Color', colors_sp(j, :), 'LineWidth', 1.5);
        
        h_legend(end+1) = h_sp;
        legend_labels{end+1} = sprintf('SP#%d d=%.4f', j, sp_results.min_distances_full(j));
    end
    
    xlabel('Species 1'); ylabel('Species 2'); zlabel('Species 3');
    title('Full Space (7D) Min Distance', 'Interpreter', 'none');
    legend(h_legend, legend_labels, 'Location', 'best', 'FontSize', 8, 'Interpreter', 'none');
    grid on; view(3);
    hold off;
    
    % 子图2: Biomass空间距离 3D相图
    subplot(2, 3, 2);
    plot3(z_traj(:, 1), z_traj(:, 2), z_traj(:, 3), 'b-', 'LineWidth', 1);
    hold on;
    
    h_legend2 = [];
    legend_labels2 = {};
    h_legend2(end+1) = plot3(NaN, NaN, NaN, 'b-', 'LineWidth', 1);
    legend_labels2{end+1} = 'Trajectory';
    
    for j = 1:n_valid_sps
        sp = sp_results.saddle_points{j};
        cp_biomass = sp_results.closest_points_biomass{j};
        
        h_sp = plot3(sp(1), sp(2), sp(3), '*', 'Color', colors_sp(j, :), ...
                     'MarkerSize', 15, 'LineWidth', 2);
        
        plot3(cp_biomass(1), cp_biomass(2), cp_biomass(3), 's', 'Color', colors_sp(j, :), ...
              'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', colors_sp(j, :));
        
        plot3([sp(1), cp_biomass(1)], [sp(2), cp_biomass(2)], [sp(3), cp_biomass(3)], ...
              '--', 'Color', colors_sp(j, :), 'LineWidth', 1.5);
        
        h_legend2(end+1) = h_sp;
        legend_labels2{end+1} = sprintf('SP#%d d=%.4f', j, sp_results.min_distances_biomass(j));
    end
    
    xlabel('Species 1'); ylabel('Species 2'); zlabel('Species 3');
    title('Biomass Space (3D) Min Distance', 'Interpreter', 'none');
    legend(h_legend2, legend_labels2, 'Location', 'best', 'FontSize', 8, 'Interpreter', 'none');
    grid on; view(3);
    hold off;
    
    % 子图3: 距离对比柱状图
    subplot(2, 3, 3);
    x_pos = 1:n_valid_sps;
    bar_width = 0.35;
    bar(x_pos - bar_width/2, sp_results.min_distances_full, bar_width, ...
        'FaceColor', [0.3, 0.6, 0.9]);
    hold on;
    bar(x_pos + bar_width/2, sp_results.min_distances_biomass, bar_width, ...
        'FaceColor', [0.9, 0.6, 0.3]);
    
    xlabel('Saddle Point Index');
    ylabel('Minimum Distance');
    title('Distance Comparison', 'Interpreter', 'none');
    legend('Full (7D)', 'Biomass (3D)', 'Location', 'best', 'Interpreter', 'none');
    xticks(1:n_valid_sps);
    xticklabels(arrayfun(@(x) sprintf('#%d', x), 1:n_valid_sps, 'UniformOutput', false));
    grid on;
    hold off;
    
    % 子图4: 全空间距离演化
    subplot(2, 3, 4);
    hold on;
    for j = 1:n_valid_sps
        sp = sp_results.saddle_points{j};
        distances_t = vecnorm(z_traj' - sp)';
        plot(t_traj, distances_t, '-', 'Color', colors_sp(j, :), 'LineWidth', 1.5);
        plot(sp_results.closest_times_full(j), sp_results.min_distances_full(j), ...
             'o', 'Color', colors_sp(j, :), 'MarkerSize', 10, 'LineWidth', 2, ...
             'MarkerFaceColor', colors_sp(j, :));
    end
    xlabel('Time');
    ylabel('Distance (Full 7D)');
    title('Full Space Distance Evolution', 'Interpreter', 'none');
    grid on;
    hold off;
    
    % 子图5: Biomass空间距离演化
    subplot(2, 3, 5);
    hold on;
    for j = 1:n_valid_sps
        sp = sp_results.saddle_points{j};
        distances_biomass_t = vecnorm(z_traj(:, 1:3)' - sp(1:3))';
        plot(t_traj, distances_biomass_t, '-', 'Color', colors_sp(j, :), 'LineWidth', 1.5);
        plot(sp_results.closest_times_biomass(j), sp_results.min_distances_biomass(j), ...
             's', 'Color', colors_sp(j, :), 'MarkerSize', 10, 'LineWidth', 2, ...
             'MarkerFaceColor', colors_sp(j, :));
    end
    xlabel('Time');
    ylabel('Distance (Biomass 3D)');
    title('Biomass Space Distance Evolution', 'Interpreter', 'none');
    grid on;
    hold off;
    
    % 子图6: 信息面板
    subplot(2, 3, 6);
    axis off;
    info_text = {
        sprintf('v123 = %.4f', v123);
        sprintf('Valid SPs: %d', n_valid_sps);
        sprintf('Biomass Var: %.4e', sp_results.biomass_var);
        '';
        'Overall Min Distances:';
        sprintf('  Full: %.6f (SP#%d)', sp_results.overall_min_dist_full, sp_results.overall_min_sp_idx_full);
        sprintf('  Biomass: %.6f (SP#%d)', sp_results.overall_min_dist_biomass, sp_results.overall_min_sp_idx_biomass);
        '';
        'Individual SP Distances:';
    };
    
    for j = 1:n_valid_sps
        info_text{end+1} = sprintf('  SP#%d: Full=%.5f, Bio=%.5f', ...
            j, sp_results.min_distances_full(j), sp_results.min_distances_biomass(j));
    end
    
    text(0.1, 0.9, info_text, 'VerticalAlignment', 'top', 'FontSize', 10, ...
         'FontName', 'Courier', 'Interpreter', 'none');
    
    sgtitle(sprintf('Parameter %d: v123 = %.4f', param_idx, v123), ...
            'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
end

function plot_three_params_comparison_biomass(results, target_params, N_spe)
    % 三参数对比图 - Biomass空间
    
    figure('Name', 'Three Parameters Comparison - Biomass Space', ...
           'Position', [100, 100, 1400, 450]);
    
    param_data = cell(3, 1);
    param_idx = 0;
    
    % 提取三个参数的数据
    for i = 1:length(results.param_results)
        if isfield(results.param_results{i}, 'has_oscillation') && ...
           results.param_results{i}.has_oscillation
            
            v123 = results.param_results{i}.v123;
            
            for j = 1:length(target_params)
                if abs(v123 - target_params(j)) < 1e-6
                    param_data{j} = results.param_results{i};
                    break;
                end
            end
        end
    end
    
    % 绘制三个子图
    for idx = 1:3
        subplot(1, 3, idx);
        
        if isempty(param_data{idx})
            text(0.5, 0.5, 'No oscillation data', 'HorizontalAlignment', 'center');
            continue;
        end
        
        sp_results = param_data{idx};
        z_traj = sp_results.trajectory_z;
        n_valid_sps = sp_results.n_valid_sps;
        colors_sp = lines(n_valid_sps);
        
        plot3(z_traj(:, 1), z_traj(:, 2), z_traj(:, 3), 'b-', 'LineWidth', 1);
        hold on;
        
        for j = 1:n_valid_sps
            sp = sp_results.saddle_points{j};
            cp_biomass = sp_results.closest_points_biomass{j};
            
            plot3(sp(1), sp(2), sp(3), '*', 'Color', colors_sp(j, :), ...
                  'MarkerSize', 15, 'LineWidth', 2);
            
            plot3(cp_biomass(1), cp_biomass(2), cp_biomass(3), 's', ...
                  'Color', colors_sp(j, :), 'MarkerSize', 10, 'LineWidth', 2, ...
                  'MarkerFaceColor', colors_sp(j, :));
            
            plot3([sp(1), cp_biomass(1)], [sp(2), cp_biomass(2)], [sp(3), cp_biomass(3)], ...
                  '--', 'Color', colors_sp(j, :), 'LineWidth', 1.5);
        end
        
        xlabel('Species 1'); ylabel('Species 2'); zlabel('Species 3');
        title(sprintf('v123 = %.4f, Min Dist = %.5f', ...
              target_params(idx), sp_results.overall_min_dist_biomass), ...
              'Interpreter', 'none');
        grid on; view(3);
        hold off;
    end
    
    sgtitle('Biomass Space (3D) - Three Parameters Comparison', ...
            'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
end

function plot_three_params_comparison_full(results, target_params, N_spe)
    % 三参数对比图 - 全空间
    
    figure('Name', 'Three Parameters Comparison - Full Space', ...
           'Position', [150, 150, 1400, 450]);
    
    param_data = cell(3, 1);
    
    % 提取三个参数的数据
    for i = 1:length(results.param_results)
        if isfield(results.param_results{i}, 'has_oscillation') && ...
           results.param_results{i}.has_oscillation
            
            v123 = results.param_results{i}.v123;
            
            for j = 1:length(target_params)
                if abs(v123 - target_params(j)) < 1e-6
                    param_data{j} = results.param_results{i};
                    break;
                end
            end
        end
    end
    
    % 绘制三个子图
    for idx = 1:3
        subplot(3, 1, idx);
        
        if isempty(param_data{idx})
            text(0.5, 0.5, 'No oscillation data', 'HorizontalAlignment', 'center');
            continue;
        end
        
        sp_results = param_data{idx};
        z_traj = sp_results.trajectory_z;
        n_valid_sps = sp_results.n_valid_sps;
        colors_sp = lines(n_valid_sps);
        
        plot3(z_traj(:, 1), z_traj(:, 2), z_traj(:, 3), 'b-', 'LineWidth', 1);
        hold on;
        
        for j = 1:n_valid_sps
            sp = sp_results.saddle_points{j};
            cp_full = sp_results.closest_points_full{j};
            
            plot3(sp(1), sp(2), sp(3), '*', 'Color', colors_sp(j, :), ...
                  'MarkerSize', 15, 'LineWidth', 2);
            
            plot3(cp_full(1), cp_full(2), cp_full(3), 'o', ...
                  'Color', colors_sp(j, :), 'MarkerSize', 10, 'LineWidth', 2, ...
                  'MarkerFaceColor', colors_sp(j, :));
            
            plot3([sp(1), cp_full(1)], [sp(2), cp_full(2)], [sp(3), cp_full(3)], ...
                  '--', 'Color', colors_sp(j, :), 'LineWidth', 1.5);
        end
        
        xlabel('[M_1]'); ylabel('[M_2]'); zlabel('[M_3]');
        title(sprintf('v123 = %.4f, Min Dist = %.5f', ...
              target_params(idx), sp_results.overall_min_dist_full), ...
              'Interpreter', 'none');
        grid on; view(3);
        hold off;
    end
    
    sgtitle('Full Space (7D) - Three Parameters Comparison', ...
            'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
end

function plot_barplot_comparison(results, target_params)
    % 柱状图对比 - 分为Biomass和Full两张图
    
    % 提取三个参数的数据
    param_data = cell(3, 1);
    for i = 1:length(results.param_results)
        if isfield(results.param_results{i}, 'has_oscillation') && ...
           results.param_results{i}.has_oscillation
            
            v123 = results.param_results{i}.v123;
            
            for j = 1:length(target_params)
                if abs(v123 - target_params(j)) < 1e-6
                    param_data{j} = results.param_results{i};
                    break;
                end
            end
        end
    end
    
    % 图1: Biomass空间柱状图
    figure('Name', 'Barplot Comparison - Biomass Space', 'Position', [200, 200, 1000, 400]);
    
    max_sps = max(cellfun(@(x) x.n_valid_sps, param_data(~cellfun(@isempty, param_data))));
    biomass_matrix = nan(3, max_sps);
    
    for i = 1:3
        if ~isempty(param_data{i})
            n_sps = param_data{i}.n_valid_sps;
            biomass_matrix(i, 1:n_sps) = param_data{i}.min_distances_biomass';
        end
    end
    
    bar(biomass_matrix');
    xlabel('Saddle Point Index');
    ylabel('Minimum Distance (Biomass 3D)');
    title('Biomass Space Distance Comparison Across Parameters', 'Interpreter', 'none');
    legend(arrayfun(@(x) sprintf('v123=%.4f', x), target_params, 'UniformOutput', false), ...
           'Location', 'best', 'Interpreter', 'none');
    grid on;
    
    % 图2: 全空间柱状图
    figure('Name', 'Barplot Comparison - Full Space', 'Position', [250, 250, 1000, 400]);
    
    full_matrix = nan(3, max_sps);
    
    for i = 1:3
        if ~isempty(param_data{i})
            n_sps = param_data{i}.n_valid_sps;
            full_matrix(i, 1:n_sps) = param_data{i}.min_distances_full';
        end
    end
    
    bar(full_matrix');
    xlabel('Saddle Point Index');
    ylabel('Minimum Distance (Full 7D)');
    title('Full Space Distance Comparison Across Parameters', 'Interpreter', 'none');
    legend(arrayfun(@(x) sprintf('v123=%.4f', x), target_params, 'UniformOutput', false), ...
           'Location', 'best', 'Interpreter', 'none');
    grid on;
end

function plot_summary_full(results, v123_list, n_params)
    % 绘制全空间距离汇总图
    
    figure('Name', 'Summary Full Space Distance', 'Position', [100, 100, 1000, 400]);
    
    % 提取有震荡的数据
    v123_osc = [];
    min_dist_full = [];
    
    for i = 1:n_params
        if isfield(results.param_results{i}, 'has_oscillation') && ...
           results.param_results{i}.has_oscillation
            v123_osc(end+1) = results.param_results{i}.v123;
            min_dist_full(end+1) = results.param_results{i}.overall_min_dist_full;
        end
    end
    
    if isempty(v123_osc)
        warning('没有检测到震荡！');
        return;
    end
    
    % 子图1: 线性尺度
    subplot(1, 2, 1);
    plot(v123_osc, min_dist_full, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('v123', 'FontSize', 12, 'Interpreter', 'none');
    ylabel('Overall Min Distance (Full 7D)', 'FontSize', 12, 'Interpreter', 'none');
    title('Full Space Distance - Linear Scale', 'FontSize', 13, 'Interpreter', 'none');
    grid on;
    
    % 子图2: 对数尺度
    subplot(1, 2, 2);
    semilogy(v123_osc, min_dist_full, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('v123', 'FontSize', 12, 'Interpreter', 'none');
    ylabel('Overall Min Distance (Full 7D, log)', 'FontSize', 12, 'Interpreter', 'none');
    title('Full Space Distance - Logarithmic Scale', 'FontSize', 13, 'Interpreter', 'none');
    grid on;
    
    sgtitle(sprintf('Full Space (7D): %d Oscillating Cases out of %d Parameters', ...
                    length(v123_osc), n_params), ...
            'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
end

function plot_summary_biomass(results, v123_list, n_params)
    % 绘制Biomass空间距离汇总图
    
    figure('Name', 'Summary Biomass Space Distance', 'Position', [150, 150, 1000, 400]);
    
    % 提取有震荡的数据
    v123_osc = [];
    min_dist_biomass = [];
    
    for i = 1:n_params
        if isfield(results.param_results{i}, 'has_oscillation') && ...
           results.param_results{i}.has_oscillation
            v123_osc(end+1) = results.param_results{i}.v123;
            min_dist_biomass(end+1) = results.param_results{i}.overall_min_dist_biomass;
        end
    end
    
    if isempty(v123_osc)
        warning('没有检测到震荡！');
        return;
    end
    
    % 子图1: 线性尺度
    subplot(1, 2, 1);
    plot(v123_osc, min_dist_biomass, 'bs-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('v123', 'FontSize', 12, 'Interpreter', 'none');
    ylabel('Overall Min Distance (Biomass 3D)', 'FontSize', 12, 'Interpreter', 'none');
    title('Biomass Space Distance - Linear Scale', 'FontSize', 13, 'Interpreter', 'none');
    grid on;
    
    % 子图2: 对数尺度
    subplot(1, 2, 2);
    semilogy(v123_osc, min_dist_biomass, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('v123', 'FontSize', 12, 'Interpreter', 'none');
    ylabel('Overall Min Distance (Biomass 3D, log)', 'FontSize', 12, 'Interpreter', 'none');
    title('Biomass Space Distance - Logarithmic Scale', 'FontSize', 13, 'Interpreter', 'none');
    grid on;
    
    sgtitle(sprintf('Biomass Space (3D): %d Oscillating Cases out of %d Parameters', ...
                    length(v123_osc), n_params), ...
            'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
end

function print_statistics(results, v123_list, n_params)
    % 打印统计信息
    
    fprintf('\n========== 统计信息 ==========\n');
    
    n_oscillating = 0;
    all_min_dists_full = [];
    all_min_dists_biomass = [];
    
    for i = 1:n_params
        if isfield(results.param_results{i}, 'has_oscillation') && ...
           results.param_results{i}.has_oscillation
            n_oscillating = n_oscillating + 1;
            all_min_dists_full(end+1) = results.param_results{i}.overall_min_dist_full;
            all_min_dists_biomass(end+1) = results.param_results{i}.overall_min_dist_biomass;
        end
    end
    
    fprintf('总参数点数: %d\n', n_params);
    fprintf('有震荡的参数点: %d (%.1f%%)\n', n_oscillating, 100*n_oscillating/n_params);
    fprintf('\n');
    
    if ~isempty(all_min_dists_full)
        fprintf('全空间(7D)最短距离统计:\n');
        fprintf('  平均值: %.6f\n', mean(all_min_dists_full));
        fprintf('  中位数: %.6f\n', median(all_min_dists_full));
        fprintf('  最小值: %.6f\n', min(all_min_dists_full));
        fprintf('  最大值: %.6f\n', max(all_min_dists_full));
        fprintf('  标准差: %.6f\n', std(all_min_dists_full));
        fprintf('\n');
        
        fprintf('Biomass空间(3D)最短距离统计:\n');
        fprintf('  平均值: %.6f\n', mean(all_min_dists_biomass));
        fprintf('  中位数: %.6f\n', median(all_min_dists_biomass));
        fprintf('  最小值: %.6f\n', min(all_min_dists_biomass));
        fprintf('  最大值: %.6f\n', max(all_min_dists_biomass));
        fprintf('  标准差: %.6f\n', std(all_min_dists_biomass));
    end
    
    fprintf('==============================\n');
end