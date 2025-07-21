function [y,D] = worst_p(k,x,R)

n = length(x);
disp(['你的3维标称向量为'])
disp(x)
disp(['TV约束距离R为 ',num2str(R)]);
disp(['你的3维系数向量为'])
disp(k)
y = zeros(3,1);

if length(k) ~= length(unique(k))
    error('向量中存在重复元素！');
end

[~,idx] = sort(k);%排序从小到大
num_3 = idx(3);
num_2 = idx(2);
num_1 = idx(1);
R_max = 2*(1-x(num_3));

alpha = min(R,R_max);
y(num_3) = x(num_3) + alpha/2;
y(num_1) = posi(x(num_1) - alpha/2);
y(num_2) = posi(x(num_2) - posi(alpha/2 - x(num_1)));

D = k'*y;
end

function y1 = posi(x1)
y1 = max(x1,0);
end


