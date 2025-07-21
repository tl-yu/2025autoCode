clc;clear
close all;

%% ��ʼ������
P_0 = [0.3 0.3 0.4;
    0.2 0.3 0.5;
    0.3 0.1 0.6];

p_0 = [0.2 0.5 0.3];%����ķֲ�
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
r = zeros(3,101); %ÿk��Ϊ��Ӧ��r(k)
u = zeros(6,100); %1��2��Ϊһģ̬�µĿ��ƣ�3��4��Ϊ��ģ̬��5��6Ϊ3ģ̬
u_pra = zeros(2,100);
x = zeros(3,101);
mode = zeros(1,101);
N = 100; %һ����101��ģ̬��ʩ��100�ο���
n = 3; %����ģ̬
worst = zeros(9,100); %100���ģ̬
u_record = zeros(1,100);
mode_record = zeros(1,3);
J = 0;

%% ��ȷ����MTPM
% ��֪�� R1, R2, R3
% ��֪�� R1, R2, R3
R1 = 0.4;
R2 = 0.8;
R3 = 0.8;
R = [R1, R2, R3];

% ʱ��
K = 100;

% ��ʼ���洢 P(k) ���еľ��� PP
PP = zeros(9, K + 1);

% ��ʼ����������Ԫ����Ϊ�ڶ�������
PP = repmat([0.5; 0.6; 0.5; 0.2; 0.3; 0.3; 0.3; 0.1; 0.2], 1, 101);

% �滻��һ��Ϊ��һ������
PP(:, 1) = [0.3; 0.2; 0.3; 0.3; 0.3; 0.1; 0.4; 0.5; 0.6];


% PP(:, 1) = P_0(:);
% 
% for k = 2:(K + 1)
%     P_k = zeros(3, 3);
%     for i = 1:3
%         % �������� TV ����Լ����ת�Ƹ�������
%         while true
%             % �����������
%             rand_vec = rand(1, 3);
%             rand_vec = rand_vec / sum(rand_vec); % ��һ��
%             % ���� TV ����
%             tv_distance = sum(abs(rand_vec - P_0(i, :)));
%             if tv_distance <= R(i)
%                 P_k(i, :) = rand_vec;
%                 break;
%             end
%         end
%     end
%     % �� P_k ����չ���洢�� PP ��
%     PP(:, k) = P_k(:);
% end

% ��ʾ���ɵ� P(k) ����
% for k = 1:(K + 1)
%     fprintf('P(%d):\n', k - 1);
%     disp(reshape(PP(:, k), 3, 3));
% end

%% ���ƹ���
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

%% ˳�ƹ���
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
        u_record(k) = 1; %��ʾ�ڵ�ǰʱ��u*û�����ڿ������ڣ�ѡ���Ե��ֵ��
    end
    J = J + q(current_mode) + w1*(u_pra(1,k))^2 + w2*(u_pra(2,k))^2;
    P_ctrl = reshape(PP(:, k+1),3,3) + u_pra(1,k)*B1' + u_pra(2,k)*B2';
    prob_vector = P_ctrl(current_mode,:);
    current_mode = choose_mode(prob_vector);
    mode(k+1) = current_mode;
end
J = J + q(current_mode)
bins = [0.5 1.5 2.5 3.5]; % �趨����߽�
mode_record = histcounts(mode(2:end), bins)
plotmode(mode);
shuzhi = sum(u_record)        
    

 
%% ������ͼ


% ����ʱ��������0-99����100���㣩
time = 0:99;

% ����ͼ�δ���
figure('Color', 'white', 'Position', [100, 100, 800, 600]);

% ���Ƶ�һ�����ݣ���ɫ��
subplot(2, 1, 1);
h1 = stairs(time, u_pra(1, :), 'Color', [0.2980 0.6863 0.3137], 'LineWidth', 1.5);
set(h1, 'LineStyle', '-', 'Marker', 'none');

% ���õ�һ��ͼ����
ax1 = gca;
ax1.XGrid = 'on';
ax1.YGrid = 'on';
ax1.GridLineStyle = '--';
ax1.XColor = [0 0 0];
ax1.YColor = [0 0 0];
ax1.FontSize = 12;
ax1.Box = 'off';

% �����������ǩ
xlabel('Time(k)', 'FontWeight', 'bold', 'Interpreter', 'latex');

% ��Ӵ��߿��ͼ�������Ͻǣ�
leg1 = legend('$u_{1k}$', 'Interpreter', 'latex', 'Location', 'northeast');
set(leg1, 'Box', 'on', 'EdgeColor', 'k', 'FontSize', 12, 'TextColor', 'k');

% ���������᷶Χ�Ϳ̶�
xlim([0 99]);
xticks(0:10:99);

% ���Ƶڶ������ݣ���ɫ��
subplot(2, 1, 2);
h2 = stairs(time, u_pra(2, :), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
set(h2, 'LineStyle', '-', 'Marker', 'none');

% ���õڶ���ͼ����
ax2 = gca;
ax2.XGrid = 'on';
ax2.YGrid = 'on';
ax2.GridLineStyle = '--';
ax2.XColor = [0 0 0];
ax2.YColor = [0 0 0];
ax2.FontSize = 12;
ax2.Box = 'off';

% �����������ǩ
xlabel('Time(k)', 'FontWeight', 'bold', 'Interpreter', 'latex');

% ��Ӵ��߿��ͼ�������Ͻǣ�
leg2 = legend('$u_{2k}$', 'Interpreter', 'latex', 'Location', 'northeast');
set(leg2, 'Box', 'on', 'EdgeColor', 'k', 'FontSize', 12, 'TextColor', 'k');

% ���������᷶Χ�Ϳ̶�
xlim([0 99]);
xticks(0:10:99);
ylim([min(u_pra(2, :)) 0.3]); % ���õڶ�����ͼ���������Ͻ�Ϊ0.1

% ʹ��subplot��Position�����ֶ�������ͼλ�ã����ֽ�С���
pos1 = get(ax1, 'Position');
pos2 = get(ax2, 'Position');
spacing = 0.05; % ���ֽ�С���
set(ax1, 'Position', [pos1(1), pos1(2) + spacing/2, pos1(3), pos1(4) - spacing/3]);
set(ax2, 'Position', [pos2(1), pos2(2), pos2(3), pos2(4) - spacing/3]);

% ʹ����ȷ��LaTeX�﷨�����ܱ���
sgtitle('$\mathbf{u_k}$', 'FontWeight', 'bold', 'Interpreter', 'latex', 'FontSize', 18); 









% % ����ʱ��������0-99����100���㣩
% time = 0:99;
% 
% % ����ͼ�δ���
% figure('Color', 'white', 'Position', [100, 100, 800, 600]);
% 
% % ������ɫ
% green = [0.2980 0.6863 0.3137];
% orange = [0.8500 0.3250 0.0980];
% black = [0 0 0];  % ��׼����
% 
% % ���Ƶ�һ������
% subplot(2, 1, 1);
% hold on;
% 
% % ����������ɫ����
% h_green = stairs(time, u_pra(1, :), 'Color', green, 'LineWidth', 1.5);
% 
% % �ڼ�¼���ϻ��ƺ�ɫ�߶�
% for k = 1:length(time)-1
%     if u_record(k) == 1
%         h_black = stairs(time(k:k+1), u_pra(1, k:k+1), 'Color', black, 'LineWidth', 1.5);
%     end
% end
% 
% % ���õ�һ��ͼ����
% ax1 = gca;
% ax1.XGrid = 'on';
% ax1.YGrid = 'on';
% ax1.GridLineStyle = '--';
% ax1.XColor = black;
% ax1.YColor = black;
% ax1.FontSize = 12;
% ax1.Box = 'off';
% 
% % �����������ǩ
% xlabel('Time(k)', 'FontWeight', 'bold', 'Interpreter', 'latex');
% ylabel('$u_{1k}$', 'Interpreter', 'latex');
% 
% % ���������귶Χ���Ͻ�Ϊ1.2��
% ylim([min(u_pra(1, :)) 1.2]);
% 
% % ��Ӵ��߿��ͼ������һ��ͼ���Ͻǣ�����λ�ñ����ص���
% leg1 = legend('$u_{1k}$ (Analytical)', '$u_{1k}$ (Vertex)', 'Interpreter', 'latex', 'Location', 'northeast');
% set(leg1, 'Box', 'on', 'EdgeColor', black, 'FontSize', 8, 'TextColor', black, 'FontWeight', 'normal');
% set(leg1, 'Position', [0.70, 0.9, 0.2, 0.04]);  % ����λ��
% leg1.NumColumns = 2;  % ����ͼ��Ϊһ������
% 
% % ��ȡͼ������Ŀ����������ɫ
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
% % ���������᷶Χ�Ϳ̶�
% xlim([0 99]);
% xticks(0:10:99);
% 
% % ���Ƶڶ�������
% subplot(2, 1, 2);
% hold on;
% 
% % ����������ɫ����
% h_orange = stairs(time, u_pra(2, :), 'Color', orange, 'LineWidth', 1.5);
% 
% % �ڼ�¼���ϻ��ƺ�ɫ�߶�
% for k = 1:length(time)-1
%     if u_record(k) == 1
%         h_black = stairs(time(k:k+1), u_pra(2, k:k+1), 'Color', black, 'LineWidth', 1.5);
%     end
% end
% 
% % ���õڶ���ͼ����
% ax2 = gca;
% ax2.XGrid = 'on';
% ax2.YGrid = 'on';
% ax2.GridLineStyle = '--';
% ax2.XColor = black;
% ax2.YColor = black;
% ax2.FontSize = 12;
% ax2.Box = 'off';
% 
% % �����������ǩ
% xlabel('Time(k)', 'FontWeight', 'bold', 'Interpreter', 'latex');
% ylabel('$u_{2k}$', 'Interpreter', 'latex');
% 
% % ���������귶Χ���Ͻ�Ϊ0.5��
% ylim([min(u_pra(2, :)) 0.2]);
% 
% % ��Ӵ��߿��ͼ�����ڶ���ͼ���Ͻǣ�����λ�ñ����ص���
% leg2 = legend('$u_{2k}$ (Analytical)', '$u_{2k}$ (Vertex)', 'Interpreter', 'latex', 'Location', 'northeast');
% set(leg2, 'Box', 'on', 'EdgeColor', black, 'FontSize', 8, 'TextColor', black, 'FontWeight', 'normal');
% set(leg2, 'Position', [0.70, 0.40, 0.2, 0.04]);  % ����λ�ã������������ͼ���ص�
% leg2.NumColumns = 2;  % ����ͼ��Ϊһ������
% 
% % ��ȡͼ������Ŀ����������ɫ
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
% % ���������᷶Χ�Ϳ̶�
% xlim([0 99]);
% xticks(0:10:99);
% 
% % ʹ��subplot��Position�����ֶ�������ͼλ�ã����ֽ�С���
% pos1 = get(ax1, 'Position');
% pos2 = get(ax2, 'Position');
% spacing = 0.05; % ���ֽ�С���
% set(ax1, 'Position', [pos1(1), pos1(2) + spacing/2, pos1(3), pos1(4) - spacing/3]);
% set(ax2, 'Position', [pos2(1), pos2(2), pos2(3), pos2(4) - spacing/3]);
% 
% % ʹ����ȷ��LaTeX�﷨�����ܱ���
% sgtitle('$\mathbf{u_k}$', 'FontWeight', 'bold', 'Interpreter', 'latex', 'FontSize', 18);
% 
% % ������Ⱦ��Ϊpainters��ȷ����ɫ��ȷ�����
% set(gcf, 'Renderer', 'painters');








% ����ʱ��������0-99����100���㣩
time = 0:99;

% ������ɫ
green = [0.2980 0.6863 0.3137];
orange = [0.8500 0.3250 0.0980];
black = [0 0 0];  % ��׼����

%% ��һ��ͼ - u_1k
figure('Color', 'white', 'Position', [100, 100, 800, 400]);
hold on;

% ����������ɫ����
h_green = stairs(time, u_pra(1, :), 'Color', green, 'LineWidth', 1.5);

% �ڼ�¼���ϻ��ƺ�ɫ�߶�
for k = 1:length(time)-1
    if u_record(k) == 1
        h_black = stairs(time(k:k+1), u_pra(1, k:k+1), 'Color', black, 'LineWidth', 1.5);
    end
end

% ��������������
ax1 = gca;
ax1.XGrid = 'on';
ax1.YGrid = 'on';
ax1.GridLineStyle = '--';
ax1.XColor = black;
ax1.YColor = black;
ax1.FontSize = 12;
ax1.Box = 'off';

% �����������ǩ
xlabel('Time(k)', 'FontWeight', 'bold', 'Interpreter', 'latex');
%ylabel('$u_{1k}$', 'Interpreter', 'latex');

% ���������귶Χ���Ͻ�Ϊ1.2��
ylim([min(u_pra(1, :))-0.1 0.7]);

% ��Ӵ��߿��ͼ�������Ͻǣ�
leg1 = legend('$u_{k}^*\in U$', '$u_{k}^*\notin U$', 'Interpreter', 'latex', 'Location', 'northeast');
% leg1 = legend('$u_{k}^*\inU_k$ ', '$u_{k}^*\notinU_k$', 'Interpreter', 'latex', 'Location', 'northeast');
set(leg1, 'Box', 'on', 'EdgeColor', black, 'FontSize', 18, 'TextColor', black, 'FontWeight', 'normal');
leg1.NumColumns = 2;  % ����ͼ��Ϊһ������

% ��ȡͼ������Ŀ����������ɫ
items1 = leg1.Children;
if length(items1) >= 2
    if isa(items1(1), 'matlab.graphics.illustration.Line')
        items1(1).Color = green;
    end
    if isa(items1(2), 'matlab.graphics.illustration.Line')
        items1(2).Color = black;
    end
end

% ���������᷶Χ�Ϳ̶�
xlim([0 100]);
xticks(0:10:100);

% ʹ����ȷ��LaTeX�﷨���ñ���
title('$u_{1k}$', 'FontWeight', 'bold', 'Interpreter', 'latex', 'FontSize', 20);

% ������Ⱦ��Ϊpainters��ȷ����ɫ��ȷ�����
set(gcf, 'Renderer', 'painters');
set(gca,'LooseInset',[0,0,0,0]);

%% �ڶ���ͼ - u_2k
% ��ȡ��Ļ�ߴ�
screenSize = get(0, 'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);

% �������λ��
figWidth = 800;
figHeight = 400;
figLeft = (screenWidth - figWidth) / 2;
figBottom = (screenHeight - figHeight) / 2;

figure('Color', 'white', 'Position', [figLeft, figBottom, figWidth, figHeight]);
hold on;

% ����������ɫ����
h_orange = stairs(time, u_pra(2, :), 'Color', orange, 'LineWidth', 1.5);

% �ڼ�¼���ϻ��ƺ�ɫ�߶�
for k = 1:length(time)-1
    if u_record(k) == 1
        h_black = stairs(time(k:k+1), u_pra(2, k:k+1), 'Color', black, 'LineWidth', 1.5);
    end
end

% ��������������
ax2 = gca;
ax2.XGrid = 'on';
ax2.YGrid = 'on';
ax2.GridLineStyle = '--';
ax2.XColor = black;
ax2.YColor = black;
ax2.FontSize = 12;
ax2.Box = 'off';

% �����������ǩ
xlabel('Time(k)', 'FontWeight', 'bold', 'Interpreter', 'latex');
%ylabel('$u_{2k}$', 'Interpreter', 'latex');

% ���������귶Χ���Ͻ�Ϊ0.2��
ylim([min(u_pra(2, :))-0.2 0.35]);

% ��Ӵ��߿��ͼ�������Ͻǣ�
leg2 = legend('$u_{k}^*\in U$', '$u_{k}^*\notin U$', 'Interpreter', 'latex', 'Location', 'northeast');
set(leg2, 'Box', 'on', 'EdgeColor', black, 'FontSize', 18, 'TextColor', black, 'FontWeight', 'normal');
leg2.NumColumns = 2;  % ����ͼ��Ϊһ������

% ��ȡͼ������Ŀ����������ɫ
items2 = leg2.Children;
if length(items2) >= 2
    if isa(items2(1), 'matlab.graphics.illustration.Line')
        items2(1).Color = orange;
    end
    if isa(items2(2), 'matlab.graphics.illustration.Line')
        items2(2).Color = black;
    end
end

% ���������᷶Χ�Ϳ̶�
xlim([0 100]);
xticks(0:10:100);

% ʹ����ȷ��LaTeX�﷨���ñ���
title('$u_{2k}$', 'FontWeight', 'bold', 'Interpreter', 'latex', 'FontSize', 20);

% ������Ⱦ��Ϊpainters��ȷ����ɫ��ȷ�����
set(gcf, 'Renderer', 'painters');    
set(gca,'LooseInset',[0,0,0,0]);