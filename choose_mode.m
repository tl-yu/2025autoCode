function mode = choose_mode(prob_vector)
% ����һ�� [0, 1] ֮��������
    r = rand;
    
    % ����������͸�������ѡ��ģ̬
    mode = find(r <= cumsum(prob_vector), 1);
end    