function [u1,u2] = optimal_cons(P,r,i)
B1 = [-0.3 0.1 0.1;
    0.2 -0.2 -0.3;
    0.1 0.1 0.2];
B2 = [0.2 0.1 0.1;
    -0.3 -0.3 0.2;
    0.1 0.2 -0.3];
A = [2 0;0 2];
B1_ = B1';
B2_ = B2';
B1_t_flat = B1_(:);  % B1ת�ú�������չ��Ϊ������
B2_t_flat = B2_(:);  % B2ת�ú�������չ��Ϊ������
b1 = r'*B1*m2s(i);
b2 = r'*B2*m2s(i);
b = [b1;b2];

% ÿ��Ԫ�ص��½�Լ����P_ij + u1*B1'_ij + u2*B2'_ij �� 0
% ת��Ϊ��׼��ʽ��B1'_ij*u1 + B2'_ij*u2 �� -P_ij
% ����-B1'_ij*u1 - B2'_ij*u2 �� P_ij
H_lower = -[B1_t_flat, B2_t_flat];
k_lower = P(:);  % ����������ӦΪ P(:)������ -P(:)

% ÿ��Ԫ�ص��Ͻ�Լ����P_ij + u1*B1'_ij + u2*B2'_ij �� 1
% ת��Ϊ��׼��ʽ��B1'_ij*u1 + B2'_ij*u2 �� 1 - P_ij
H_upper = [B1_t_flat, B2_t_flat];
k_upper = 1 - P(:);

% �ϲ����½�Լ��
H = [H_lower; H_upper];  % 2n?��2 ����
k = [k_lower; k_upper];  % 2n?��1 ����

% ʹ��quadprog�����ι滮����
options = optimoptions('quadprog', 'Display', 'iter');
u_star = quadprog(A, b, H, k, [], [], [], [], [], options);

% �������Ž��Ӧ��ת�ƾ���
PP_star = P + u_star(1)*B1' + u_star(2)*B2';

% ��ʾ���
fprintf('���Ž� u* = [%.4f, %.4f]''\n', u_star(1), u_star(2));
fprintf('��֤ PP* �к�: ');
fprintf('%.4f ', sum(PP_star, 2)');
fprintf('\n');

% ��֤Լ������
constraint_values = H * u_star - k;
max_violation = max(constraint_values);
fprintf('���Լ��Υ��ֵ: %.4e\n', max_violation);
if max_violation <= 1e-8
    fprintf('OK ����Լ�����㣨��Χ�ڣ�\n');
else
    fprintf('No ����Լ��Υ��\n');
    % ��ʾΥ��Լ����Ԫ��λ��
    violating_indices = find(constraint_values > 1e-8);
    fprintf('Υ��Լ����Ԫ������: ');
    fprintf('%d ', violating_indices');
    fprintf('\n');
end
u1 = u_star(1);
u2 = u_star(2);

end