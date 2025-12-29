% 位置生成函数
function positions = generate_positions(N, lambda)
    % 参数说明:
    % N: 天线数量
    % lambda: 波长
    % 返回值:
    % positions: 生成的天线位置 (2 x N)

    % 最小间距
    min_spacing = lambda / 2;

    % 初始化位置矩阵
    positions = zeros(2, N);

    % 从中心点开始向外扩散
    center = [0; 0]; % 中心点
    positions(:, 1) = center; % 第一个天线在中心

    % 生成其余天线的位置
    for i = 2:N
        while true
            % 随机生成一个候选位置
            candidate = (rand(2, 1) - 0.5) * 3 * lambda;

            % 检查候选位置是否满足最小间距要求
            valid = true;
            for j = 1:i-1
                if norm(candidate - positions(:, j)) < min_spacing
                    valid = false;
                    break;
                end
            end

            % 如果满足要求，则保存候选位置
            if valid
                positions(:, i) = candidate;
                break;
            end
        end
    end
end