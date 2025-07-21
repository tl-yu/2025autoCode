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
B1_t_flat = B1_(:);  % B1转置后按列优先展开为列向量
B2_t_flat = B2_(:);  % B2转置后按列优先展开为列向量
b1 = r'*B1*m2s(i);
b2 = r'*B2*m2s(i);
b = [b1;b2];

% 每个元素的下界约束：P_ij + u1*B1'_ij + u2*B2'_ij ≥ 0
% 转化为标准形式：B1'_ij*u1 + B2'_ij*u2 ≥ -P_ij
% 即：-B1'_ij*u1 - B2'_ij*u2 ≤ P_ij
H_lower = -[B1_t_flat, B2_t_flat];
k_lower = P(:);  % 修正：这里应为 P(:)，而非 -P(:)

% 每个元素的上界约束：P_ij + u1*B1'_ij + u2*B2'_ij ≤ 1
% 转化为标准形式：B1'_ij*u1 + B2'_ij*u2 ≤ 1 - P_ij
H_upper = [B1_t_flat, B2_t_flat];
k_upper = 1 - P(:);

% 合并上下界约束
H = [H_lower; H_upper];  % 2n?×2 矩阵
k = [k_lower; k_upper];  % 2n?×1 向量

% 使用quadprog求解二次规划问题
options = optimoptions('quadprog', 'Display', 'iter');
u_star = quadprog(A, b, H, k, [], [], [], [], [], options);

% 计算最优解对应的转移矩阵
PP_star = P + u_star(1)*B1' + u_star(2)*B2';

% 显示结果
fprintf('最优解 u* = [%.4f, %.4f]''\n', u_star(1), u_star(2));
fprintf('验证 PP* 行和: ');
fprintf('%.4f ', sum(PP_star, 2)');
fprintf('\n');

% 验证约束条件
constraint_values = H * u_star - k;
max_violation = max(constraint_values);
fprintf('最大约束违反值: %.4e\n', max_violation);
if max_violation <= 1e-8
    fprintf('OK 所有约束满足（误差范围内）\n');
else
    fprintf('No 存在约束违反\n');
    % 显示违反约束的元素位置
    violating_indices = find(constraint_values > 1e-8);
    fprintf('违反约束的元素索引: ');
    fprintf('%d ', violating_indices');
    fprintf('\n');
end
u1 = u_star(1);
u2 = u_star(2);

end