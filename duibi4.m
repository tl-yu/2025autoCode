clc;clear

mode_total = zeros(200,3);
J_total = 0;

for ite = 1:200
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
    
    
%     PP(:, 1) = P_0(:);
% 
%     for k = 2:(K + 1)
%         P_k = zeros(3, 3);
%         for i = 1:3
%             % �������� TV ����Լ����ת�Ƹ�������
%             while true
%                 % �����������
%                 rand_vec = rand(1, 3);
%                 rand_vec = rand_vec / sum(rand_vec); % ��һ��
%                 % ���� TV ����
%                 tv_distance = sum(abs(rand_vec - P_0(i, :)));
%                 if tv_distance <= R(i)
%                     P_k(i, :) = rand_vec;
%                     break;
%                 end
%             end
%         end
%         % �� P_k ����չ���洢�� PP ��
%         PP(:, k) = P_k(:);
%     end

    % ��ʾ���ɵ� P(k) ����
    % for k = 1:(K + 1)
    %     fprintf('P(%d):\n', k - 1);
    %     disp(reshape(PP(:, k), 3, 3));
    % end


    %% ˳�ƹ���
    temp1 = rand;
    current_mode = find(temp1 <= cumsum(p_0), 1);
    mode(1) = current_mode;
    for k = 1:100
        u_pra(1,k) = 2*rand - 1;
        u_pra(2,k) = 2*rand - 1;
        while inU(reshape(PP(:, k+1), 3, 3),u_pra(1,k),u_pra(2,k)) ~= 1
                u_pra(1,k) = 2*rand - 1;
                u_pra(2,k) = 2*rand - 1;
        end
        disp('���ѡ��Ŀ�����u���ڿ�������')
        J = J + q(current_mode) + w1*(u_pra(1,k))^2 + w2*(u_pra(2,k))^2;
        P_ctrl = reshape(PP(:, k+1),3,3) + u_pra(1,k)*B1' + u_pra(2,k)*B2';
        prob_vector = P_ctrl(current_mode,:);
        current_mode = choose_mode(prob_vector);
        mode(k+1) = current_mode;
    end
    J = J + q(current_mode)
    bins = [0.5 1.5 2.5 3.5]; % �趨����߽�
    mode_record = histcounts(mode(2:end), bins)
    J_total = J_total + J;
    mode_total(ite,:) = mode_record;
    
end

disp('---------------------------------------------')
disp('��4������µĴ�����ģ̬�ֲ�')   
J_total
fenbu = sum(mode_total)
fenbu_p = fenbu./(100*200)   


