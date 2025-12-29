% 参数设置
K = 4; % 用户数量
M = 8; % 总路径数（包括直射路径和散射路径）
user_center = [40; 30; 0]; % 用户区域中心
radius = 10; % 用户区域半径
GBS_pos = [0; 0; 30]; % GBS 位置
LBS_pos = [0; 100; 30]; % LBS 位置
UAV_pos = [10; 20; 100]; % 无人机位置
Jammer_pos = [40; 80; 0]; % 干扰源位置
output_file = 'dataset2.mat'; % 输出文件名

% 调用函数
generate_channel_param(K, M, user_center, radius, GBS_pos, LBS_pos, UAV_pos, Jammer_pos, output_file);

% load('communication_param.mat')