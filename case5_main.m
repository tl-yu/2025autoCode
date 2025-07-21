clc;clear
close all;

%% 初始化参数
P_0 = [0.3 0.3 0.4;
    0.2 0.3 0.5;
    0.3 0.1 0.6];

p_0 = [0.2 0.5 0.3];%最初的分布
q = [7,4,1]';
w1 = 1, w2 = 1;
W = diag([w1 w2]);
f = q;
B1 = [-0.3 0.1 0.1;
    0.2 -0.2 -0.3;
    0.1 0.1 0.2];
B2 = [0.2 0.1 0.1;
    -0.3 -0.3 0.2;
    0.1 0.2 -0.3];
r = zeros(3,101); %每k列为对应的r(k)
u = zeros(6,100); %1，2行为一模态下的控制，3，4行为二模态，5，6为3模态
u_pra = zeros(2,100);
x = zeros(3,101);
mode = zeros(1,101);
N = 100; %一共有101个模态，施加100次控制
n = 3; %三个模态
worst = zeros(9,100); %100个最坏模态
u_record = zeros(1,100);
mode_record = zeros(1,3);
J = 0;

%% 不确定的MTPM
% 已知的 R1, R2, R3
% 已知的 R1, R2, R3
R1 = 0.4;
R2 = 0.8;
R3 = 0.8;
R = [R1, R2, R3];

% 时步
K = 100;

% 初始化存储 P(k) 序列的矩阵 PP
PP = zeros(9, K + 1);

% 初始化矩阵，所有元素设为第二组数据
PP = repmat([0.5; 0.6; 0.5; 0.2; 0.3; 0.3; 0.3; 0.1; 0.2], 1, 101);

% 替换第一列为第一组数据
PP(:, 1) = [0.3; 0.2; 0.3; 0.3; 0.3; 0.1; 0.4; 0.5; 0.6];


% PP(:, 1) = P_0(:);
% 
% for k = 2:(K + 1)
%     P_k = zeros(3, 3);
%     for i = 1:3
%         % 生成满足 TV 距离约束的转移概率向量
%         while true
%             % 生成随机向量
%             rand_vec = rand(1, 3);
%             rand_vec = rand_vec / sum(rand_vec); % 归一化
%             % 计算 TV 距离
%             tv_distance = sum(abs(rand_vec - P_0(i, :)));
%             if tv_distance <= R(i)
%                 P_k(i, :) = rand_vec;
%                 break;
%             end
%         end
%     end
%     % 将 P_k 按列展开存储到 PP 中
%     PP(:, k) = P_k(:);
% end

% 显示生成的 P(k) 序列
% for k = 1:(K + 1)
%     fprintf('P(%d):\n', k - 1);
%     disp(reshape(PP(:, k), 3, 3));
% end

%% 倒推过程
r(:,end) = f;
for k = 100:-1:1
    for i = 1:3
        u(2*i-1,k) = -(r(:,k+1)'*B1*m2s(i))/(2*w1);
        u(2*i,k) = -(r(:,k+1)'*B2*m2s(i))/(2*w2);
        p(:,i) = worst_p(r(:,k+1),P_0(i,:)',R(i));
        worst(3*i-2:3*i,k) = p(:,i);
        r(i,k) = q(i) + r(:,k+1)'*p(:,i) - (r(:,k+1)'*B1*m2s(i))/(4*w1) - (r(:,k+1)'*B2*m2s(i))/(4*w2);
    end
end

%% 顺推过程
temp1 = rand;
current_mode = find(temp1 <= cumsum(p_0), 1);
mode(1) = current_mode;
for k = 1:100
    if inU(reshape(PP(:, k+1), 3, 3),u(2*current_mode-1,k),u(2*current_mode,k)) == 1
        u_pra(1,k) = u(2*current_mode-1,k);
        u_pra(2,k) = u(2*current_mode,k);
    else
        [uu1,uu2] = optimal_cons(reshape(PP(:, k+1),3,3),r(:,k+1),current_mode);
        u_pra(1,k) = uu1;
        u_pra(2,k) = uu2;
        u_record(k) = 1; %表示在当前时步u*没有落在可行域内，选择边缘数值解
    end
    J = J + q(current_mode) + w1*(u_pra(1,k))^2 + w2*(u_pra(2,k))^2;
    P_ctrl = reshape(PP(:, k+1),3,3) + u_pra(1,k)*B1' + u_pra(2,k)*B2';
    prob_vector = P_ctrl(current_mode,:);
    current_mode = choose_mode(prob_vector);
    mode(k+1) = current_mode;
end
J = J + q(current_mode)
bins = [0.5 1.5 2.5 3.5]; % 设定区间边界
mode_record = histcounts(mode(2:end), bins)
plotmode(mode);
shuzhi = sum(u_record)        
    

 
%% 画控制图


% 创建时间向量（0-99，共100个点）
time = 0:99;

% 创建图形窗口
figure('Color', 'white', 'Position', [100, 100, 800, 600]);

% 绘制第一行数据（绿色）
subplot(2, 1, 1);
h1 = stairs(time, u_pra(1, :), 'Color', [0.2980 0.6863 0.3137], 'LineWidth', 1.5);
set(h1, 'LineStyle', '-', 'Marker', 'none');

% 设置第一子图属性
ax1 = gca;
ax1.XGrid = 'on';
ax1.YGrid = 'on';
ax1.GridLineStyle = '--';
ax1.XColor = [0 0 0];
ax1.YColor = [0 0 0];
ax1.FontSize = 12;
ax1.Box = 'off';

% 设置坐标轴标签
xlabel('Time(k)', 'FontWeight', 'bold', 'Interpreter', 'latex');

% 添加带边框的图例（右上角）
leg1 = legend('$u_{1k}$', 'Interpreter', 'latex', 'Location', 'northeast');
set(leg1, 'Box', 'on', 'EdgeColor', 'k', 'FontSize', 12, 'TextColor', 'k');

% 设置坐标轴范围和刻度
xlim([0 99]);
xticks(0:10:99);

% 绘制第二行数据（橙色）
subplot(2, 1, 2);
h2 = stairs(time, u_pra(2, :), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
set(h2, 'LineStyle', '-', 'Marker', 'none');

% 设置第二子图属性
ax2 = gca;
ax2.XGrid = 'on';
ax2.YGrid = 'on';
ax2.GridLineStyle = '--';
ax2.XColor = [0 0 0];
ax2.YColor = [0 0 0];
ax2.FontSize = 12;
ax2.Box = 'off';

% 设置坐标轴标签
xlabel('Time(k)', 'FontWeight', 'bold', 'Interpreter', 'latex');

% 添加带边框的图例（右上角）
leg2 = legend('$u_{2k}$', 'Interpreter', 'latex', 'Location', 'northeast');
set(leg2, 'Box', 'on', 'EdgeColor', 'k', 'FontSize', 12, 'TextColor', 'k');

% 设置坐标轴范围和刻度
xlim([0 99]);
xticks(0:10:99);
ylim([min(u_pra(2, :)) 0.3]); % 设置第二个子图的纵坐标上界为0.1

% 使用subplot的Position属性手动调整子图位置，保持较小间距
pos1 = get(ax1, 'Position');
pos2 = get(ax2, 'Position');
spacing = 0.05; % 保持较小间距
set(ax1, 'Position', [pos1(1), pos1(2) + spacing/2, pos1(3), pos1(4) - spacing/3]);
set(ax2, 'Position', [pos2(1), pos2(2), pos2(3), pos2(4) - spacing/3]);

% 使用正确的LaTeX语法设置总标题
sgtitle('$\mathbf{u_k}$', 'FontWeight', 'bold', 'Interpreter', 'latex', 'FontSize', 18); 









% % 创建时间向量（0-99，共100个点）
% time = 0:99;
% 
% % 创建图形窗口
% figure('Color', 'white', 'Position', [100, 100, 800, 600]);
% 
% % 定义颜色
% green = [0.2980 0.6863 0.3137];
% orange = [0.8500 0.3250 0.0980];
% black = [0 0 0];  % 标准纯黑
% 
% % 绘制第一行数据
% subplot(2, 1, 1);
% hold on;
% 
% % 绘制完整绿色线条
% h_green = stairs(time, u_pra(1, :), 'Color', green, 'LineWidth', 1.5);
% 
% % 在记录点上绘制黑色线段
% for k = 1:length(time)-1
%     if u_record(k) == 1
%         h_black = stairs(time(k:k+1), u_pra(1, k:k+1), 'Color', black, 'LineWidth', 1.5);
%     end
% end
% 
% % 设置第一子图属性
% ax1 = gca;
% ax1.XGrid = 'on';
% ax1.YGrid = 'on';
% ax1.GridLineStyle = '--';
% ax1.XColor = black;
% ax1.YColor = black;
% ax1.FontSize = 12;
% ax1.Box = 'off';
% 
% % 设置坐标轴标签
% xlabel('Time(k)', 'FontWeight', 'bold', 'Interpreter', 'latex');
% ylabel('$u_{1k}$', 'Interpreter', 'latex');
% 
% % 设置纵坐标范围（上界为1.2）
% ylim([min(u_pra(1, :)) 1.2]);
% 
% % 添加带边框的图例（第一子图右上角，调整位置避免重叠）
% leg1 = legend('$u_{1k}$ (Analytical)', '$u_{1k}$ (Vertex)', 'Interpreter', 'latex', 'Location', 'northeast');
% set(leg1, 'Box', 'on', 'EdgeColor', black, 'FontSize', 8, 'TextColor', black, 'FontWeight', 'normal');
% set(leg1, 'Position', [0.70, 0.9, 0.2, 0.04]);  % 调整位置
% leg1.NumColumns = 2;  % 设置图例为一行两列
% 
% % 获取图例的条目对象并设置颜色
% items1 = leg1.Children;
% if length(items1) >= 2
%     if isa(items1(1), 'matlab.graphics.illustration.Line')
%         items1(1).Color = green;
%     end
%     if isa(items1(2), 'matlab.graphics.illustration.Line')
%         items1(2).Color = black;
%     end
% end
% 
% % 设置坐标轴范围和刻度
% xlim([0 99]);
% xticks(0:10:99);
% 
% % 绘制第二行数据
% subplot(2, 1, 2);
% hold on;
% 
% % 绘制完整橙色线条
% h_orange = stairs(time, u_pra(2, :), 'Color', orange, 'LineWidth', 1.5);
% 
% % 在记录点上绘制黑色线段
% for k = 1:length(time)-1
%     if u_record(k) == 1
%         h_black = stairs(time(k:k+1), u_pra(2, k:k+1), 'Color', black, 'LineWidth', 1.5);
%     end
% end
% 
% % 设置第二子图属性
% ax2 = gca;
% ax2.XGrid = 'on';
% ax2.YGrid = 'on';
% ax2.GridLineStyle = '--';
% ax2.XColor = black;
% ax2.YColor = black;
% ax2.FontSize = 12;
% ax2.Box = 'off';
% 
% % 设置坐标轴标签
% xlabel('Time(k)', 'FontWeight', 'bold', 'Interpreter', 'latex');
% ylabel('$u_{2k}$', 'Interpreter', 'latex');
% 
% % 设置纵坐标范围（上界为0.5）
% ylim([min(u_pra(2, :)) 0.2]);
% 
% % 添加带边框的图例（第二子图右上角，调整位置避免重叠）
% leg2 = legend('$u_{2k}$ (Analytical)', '$u_{2k}$ (Vertex)', 'Interpreter', 'latex', 'Location', 'northeast');
% set(leg2, 'Box', 'on', 'EdgeColor', black, 'FontSize', 8, 'TextColor', black, 'FontWeight', 'normal');
% set(leg2, 'Position', [0.70, 0.40, 0.2, 0.04]);  % 调整位置，避免与上面的图例重叠
% leg2.NumColumns = 2;  % 设置图例为一行两列
% 
% % 获取图例的条目对象并设置颜色
% items2 = leg2.Children;
% if length(items2) >= 2
%     if isa(items2(1), 'matlab.graphics.illustration.Line')
%         items2(1).Color = orange;
%     end
%     if isa(items2(2), 'matlab.graphics.illustration.Line')
%         items2(2).Color = black;
%     end
% end
% 
% % 设置坐标轴范围和刻度
% xlim([0 99]);
% xticks(0:10:99);
% 
% % 使用subplot的Position属性手动调整子图位置，保持较小间距
% pos1 = get(ax1, 'Position');
% pos2 = get(ax2, 'Position');
% spacing = 0.05; % 保持较小间距
% set(ax1, 'Position', [pos1(1), pos1(2) + spacing/2, pos1(3), pos1(4) - spacing/3]);
% set(ax2, 'Position', [pos2(1), pos2(2), pos2(3), pos2(4) - spacing/3]);
% 
% % 使用正确的LaTeX语法设置总标题
% sgtitle('$\mathbf{u_k}$', 'FontWeight', 'bold', 'Interpreter', 'latex', 'FontSize', 18);
% 
% % 设置渲染器为painters（确保颜色精确输出）
% set(gcf, 'Renderer', 'painters');








% 创建时间向量（0-99，共100个点）
time = 0:99;

% 定义颜色
green = [0.2980 0.6863 0.3137];
orange = [0.8500 0.3250 0.0980];
black = [0 0 0];  % 标准纯黑

%% 第一张图 - u_1k
figure('Color', 'white', 'Position', [100, 100, 800, 400]);
hold on;

% 绘制完整绿色线条
h_green = stairs(time, u_pra(1, :), 'Color', green, 'LineWidth', 1.5);

% 在记录点上绘制黑色线段
for k = 1:length(time)-1
    if u_record(k) == 1
        h_black = stairs(time(k:k+1), u_pra(1, k:k+1), 'Color', black, 'LineWidth', 1.5);
    end
end

% 设置坐标轴属性
ax1 = gca;
ax1.XGrid = 'on';
ax1.YGrid = 'on';
ax1.GridLineStyle = '--';
ax1.XColor = black;
ax1.YColor = black;
ax1.FontSize = 12;
ax1.Box = 'off';

% 设置坐标轴标签
xlabel('Time(k)', 'FontWeight', 'bold', 'Interpreter', 'latex');
%ylabel('$u_{1k}$', 'Interpreter', 'latex');

% 设置纵坐标范围（上界为1.2）
ylim([min(u_pra(1, :))-0.1 0.7]);

% 添加带边框的图例（右上角）
leg1 = legend('$u_{k}^*\in U$', '$u_{k}^*\notin U$', 'Interpreter', 'latex', 'Location', 'northeast');
% leg1 = legend('$u_{k}^*\inU_k$ ', '$u_{k}^*\notinU_k$', 'Interpreter', 'latex', 'Location', 'northeast');
set(leg1, 'Box', 'on', 'EdgeColor', black, 'FontSize', 18, 'TextColor', black, 'FontWeight', 'normal');
leg1.NumColumns = 2;  % 设置图例为一行两列

% 获取图例的条目对象并设置颜色
items1 = leg1.Children;
if length(items1) >= 2
    if isa(items1(1), 'matlab.graphics.illustration.Line')
        items1(1).Color = green;
    end
    if isa(items1(2), 'matlab.graphics.illustration.Line')
        items1(2).Color = black;
    end
end

% 设置坐标轴范围和刻度
xlim([0 100]);
xticks(0:10:100);

% 使用正确的LaTeX语法设置标题
title('$u_{1k}$', 'FontWeight', 'bold', 'Interpreter', 'latex', 'FontSize', 20);

% 设置渲染器为painters（确保颜色精确输出）
set(gcf, 'Renderer', 'painters');
set(gca,'LooseInset',[0,0,0,0]);

%% 第二张图 - u_2k
% 获取屏幕尺寸
screenSize = get(0, 'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);

% 计算居中位置
figWidth = 800;
figHeight = 400;
figLeft = (screenWidth - figWidth) / 2;
figBottom = (screenHeight - figHeight) / 2;

figure('Color', 'white', 'Position', [figLeft, figBottom, figWidth, figHeight]);
hold on;

% 绘制完整橙色线条
h_orange = stairs(time, u_pra(2, :), 'Color', orange, 'LineWidth', 1.5);

% 在记录点上绘制黑色线段
for k = 1:length(time)-1
    if u_record(k) == 1
        h_black = stairs(time(k:k+1), u_pra(2, k:k+1), 'Color', black, 'LineWidth', 1.5);
    end
end

% 设置坐标轴属性
ax2 = gca;
ax2.XGrid = 'on';
ax2.YGrid = 'on';
ax2.GridLineStyle = '--';
ax2.XColor = black;
ax2.YColor = black;
ax2.FontSize = 12;
ax2.Box = 'off';

% 设置坐标轴标签
xlabel('Time(k)', 'FontWeight', 'bold', 'Interpreter', 'latex');
%ylabel('$u_{2k}$', 'Interpreter', 'latex');

% 设置纵坐标范围（上界为0.2）
ylim([min(u_pra(2, :))-0.2 0.35]);

% 添加带边框的图例（右上角）
leg2 = legend('$u_{k}^*\in U$', '$u_{k}^*\notin U$', 'Interpreter', 'latex', 'Location', 'northeast');
set(leg2, 'Box', 'on', 'EdgeColor', black, 'FontSize', 18, 'TextColor', black, 'FontWeight', 'normal');
leg2.NumColumns = 2;  % 设置图例为一行两列

% 获取图例的条目对象并设置颜色
items2 = leg2.Children;
if length(items2) >= 2
    if isa(items2(1), 'matlab.graphics.illustration.Line')
        items2(1).Color = orange;
    end
    if isa(items2(2), 'matlab.graphics.illustration.Line')
        items2(2).Color = black;
    end
end

% 设置坐标轴范围和刻度
xlim([0 100]);
xticks(0:10:100);

% 使用正确的LaTeX语法设置标题
title('$u_{2k}$', 'FontWeight', 'bold', 'Interpreter', 'latex', 'FontSize', 20);

% 设置渲染器为painters（确保颜色精确输出）
set(gcf, 'Renderer', 'painters');    
set(gca,'LooseInset',[0,0,0,0]);