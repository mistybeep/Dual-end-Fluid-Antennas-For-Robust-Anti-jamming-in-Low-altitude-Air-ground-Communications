function movable_coordinates = generate_upa_coordinates(reg_length, num_movable_antennas, wave_length)
    % 生成正方形区域中天线的最优摆放位置，使天线最大化覆盖正方形区域。
    % 输入:
    %   reg_length: 正方形区域长度的倍数
    %   num_movable_antennas: 天线的数量（即圆的数量）
    %   wave_length: 波长，用于确定半波长
    % 输出:
    %   movable_coordinates: 天线的坐标列表

    % 参数设置
    lambda_half = wave_length / 2;  % 半波长
    max_boundary = reg_length * lambda_half;  % 正方形边界长度
    TOL = 1e-10;  % 提高数值稳定性

    % 初始猜测：均匀分布圆心，半径取边界长度的合理初始值
    x0 = zeros(2*num_movable_antennas + 1, 1);
    grid_size = ceil(sqrt(num_movable_antennas));  % 按网格初步分布
    spacing = 2*max_boundary / (grid_size + 1);
    idx = 1;
    for i = 1:grid_size
        for j = 1:grid_size
            if idx <= num_movable_antennas
                x0(2*idx-1) = -max_boundary + j*spacing;  % x 坐标
                x0(2*idx) = -max_boundary + i*spacing;    % y 坐标
                idx = idx + 1;
            end
        end
    end
    x0(end) = max_boundary / (2*sqrt(num_movable_antennas));  % 初始半径猜测

    % 目标函数：最大化半径（最小化负半径）
    obj_fun = @(z) -z(end);

    % 非线性约束
    nonlcon = @(z) nonlinear_constraints(z, num_movable_antennas, max_boundary, TOL);

    % 上下界
    lb = [-max_boundary*ones(2*num_movable_antennas, 1); 0];  % x, y 下界和 r >= 0
    ub = [max_boundary*ones(2*num_movable_antennas, 1); max_boundary];  % r 有上界

    % 优化选项
    options = optimoptions('fmincon', ...
        'Algorithm', 'interior-point', ...
        'MaxIterations', 1000, ...
        'TolFun', 1e-8, ...
        'TolCon', 1e-8, ...
        'Display', 'off');  % 关闭显示
    
    % 求解优化问题
    [sol, fval, exitflag] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, nonlcon, options);

    % 检查求解状态
    % if exitflag <= 0
    %     warning('优化未成功收敛，exitflag = %d', exitflag);
    % end

    % 提取结果
    movable_coordinates = zeros(2, num_movable_antennas);
    for i = 1:num_movable_antennas
        movable_coordinates(1,i) = real(sol(2*i - 1));
        movable_coordinates(2,i) = real(sol(2*i));
    end
    % movable_coordinates = sol(1:end-1);  % x 和 y 坐标
    r_opt = sol(end);  % 最优半径

    % 输出结果验证
    % fprintf('最优半径: %.4f\n', r_opt);
    % fprintf('天线坐标:\n');
    % for i = 1:num_movable_antennas
    %     fprintf('天线 %d: (%.4f, %.4f)\n', i, movable_coordinates(2*i-1), movable_coordinates(2*i));
    % end

    % 绘图验证
    % figure;
    % hold on;
    % axis equal;
    % theta = linspace(0, 2*pi, 100);
    % for i = 1:num_movable_antennas
    %     x_center = movable_coordinates(2*i-1);
    %     y_center = movable_coordinates(2*i);
    %     x_circle = x_center + r_opt*cos(theta);
    %     y_circle = y_center + r_opt*sin(theta);
    %     plot(x_circle, y_circle, 'b-', 'LineWidth', 1.5);
    %     plot(x_center, y_center, 'ro');  % 标记圆心
    % end
    % xlim([-max_boundary max_boundary]);
    % ylim([-max_boundary max_boundary]);
    % rectangle('Position', [-max_boundary -max_boundary 2*max_boundary 2*max_boundary], ...
    %           'LineStyle', '--', 'LineWidth', 1);
    % title(sprintf('Optimal Antenna Placement (r = %.4f)', r_opt));
    % grid on;
    % hold off;

end

function [c, ceq] = nonlinear_constraints(z, num_movable_antennas, max_boundary, TOL)
    % 非线性约束函数
    x = z(1:2:end-1);  % x 坐标
    y = z(2:2:end-1);  % y 坐标
    r = z(end);  % 半径

    % 1. 每个圆必须在正方形内
    c_boundary = [
        x + r - max_boundary;    % 右侧边界
        y + r - max_boundary;    % 上侧边界
        -x + r - max_boundary;   % 左侧边界
        -y + r - max_boundary    % 下侧边界
    ];

    % 2. 每两个圆之间不能重叠
    c_overlap = [];
    for i = 1:num_movable_antennas
        for j = i+1:num_movable_antennas
            dist = sqrt((x(i) - x(j))^2 + (y(i) - y(j))^2 + TOL);  % 避免零距离
            c_overlap = [c_overlap; 2*r - dist];
        end
    end

    % 所有不等式约束 (c <= 0)
    c = [c_boundary; c_overlap];
    ceq = [];  % 无等式约束
end