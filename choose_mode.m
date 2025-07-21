function mode = choose_mode(prob_vector)
% 生成一个 [0, 1] 之间的随机数
    r = rand;
    
    % 根据随机数和概率向量选择模态
    mode = find(r <= cumsum(prob_vector), 1);
end    