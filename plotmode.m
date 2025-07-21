function plotmode(mode_vector)
    % ��֤�����Ƿ�Ϊ����
    if ~isvector(mode_vector)
        error('�������������');
    end
    
    % ��׼��ʱ�䲽����0-100
    time_steps = length(mode_vector);
    time = linspace(0, 100, time_steps);
    
    % ����ͼ�δ���
    figure('Color', 'white');
    
    % ���ƽ���״ͼ����������
    h = stairs(time, mode_vector, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
    set(h, 'LineStyle', '-', 'Marker', 'none');
    
    % ��������������
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridLineStyle = '--';
    ax.XColor = [0 0 0];
    ax.YColor = [0 0 0];
    ax.FontSize = 12;
    ax.Box = 'off'; % �Ƴ��߿�
    
    % �����������ǩ�ͱ���
    xlabel('Time(k)', 'FontWeight', 'bold','Interpreter', 'latex');
    ylabel('Mode $\theta_k$', 'Interpreter', 'latex', 'FontWeight', 'bold');
%     title('Mode Over Time', 'FontWeight', 'bold','Interpreter', 'latex');
    
    % ����x�᷶Χ��0��100���޶���հ�
    xlim([0 100]);
    
    % ���ÿ̶�
    xticks(0:10:100);
    ylim([0 4]);        % �̶�Y����ʾ��ΧΪ0-4
    yticks(0:1:4);      % ���ÿ̶ȴ�0��4������Ϊ1
    set(gca,'LooseInset',[0,0,0,0]);
end