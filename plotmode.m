function plotmode(mode_vector)
    % 验证输入是否为向量
    if ~isvector(mode_vector)
        error('输入必须是向量');
    end
    
    % 标准化时间步长到0-100
    time_steps = length(mode_vector);
    time = linspace(0, 100, time_steps);
    
    % 创建图形窗口
    figure('Color', 'white');
    
    % 绘制阶梯状图（右连续）
    h = stairs(time, mode_vector, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
    set(h, 'LineStyle', '-', 'Marker', 'none');
    
    % 设置坐标轴属性
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridLineStyle = '--';
    ax.XColor = [0 0 0];
    ax.YColor = [0 0 0];
    ax.FontSize = 12;
    ax.Box = 'off'; % 移除边框
    
    % 设置坐标轴标签和标题
    xlabel('Time(k)', 'FontWeight', 'bold','Interpreter', 'latex');
    ylabel('Mode $\theta_k$', 'Interpreter', 'latex', 'FontWeight', 'bold');
%     title('Mode Over Time', 'FontWeight', 'bold','Interpreter', 'latex');
    
    % 设置x轴范围从0到100，无多余空白
    xlim([0 100]);
    
    % 设置刻度
    xticks(0:10:100);
    ylim([0 4]);        % 固定Y轴显示范围为0-4
    yticks(0:1:4);      % 设置刻度从0到4，步长为1
    set(gca,'LooseInset',[0,0,0,0]);
end