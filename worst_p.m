function [y,D] = worst_p(k,x,R)

n = length(x);
disp(['���3ά�������Ϊ'])
disp(x)
disp(['TVԼ������RΪ ',num2str(R)]);
disp(['���3άϵ������Ϊ'])
disp(k)
y = zeros(3,1);

if length(k) ~= length(unique(k))
    error('�����д����ظ�Ԫ�أ�');
end

[~,idx] = sort(k);%�����С����
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


