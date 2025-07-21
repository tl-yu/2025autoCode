function [y,D,n,i,r] = worst_rate(x,k,R)

%% 预设置
n = length(x);
disp(['你的 ',num2str(n),' 维标称向量为'])
disp(x)
i = find(x < 0);
disp(['当前状态向量i为 ',num2str(i)]);
disp(['TV约束距离R(i)为 ',num2str(R)]);
disp(['你的 ',num2str(n),' 维系数向量为'])
disp(k)

%% 划分状态集 
% 确保输入是非负向量
if any(k < 0)
    error('输入向量 k 必须是非负的');
end

% 提取唯一值并排序
[~, ~, unique_indices] = unique(k); % 提取唯一值索引
sorted_ranks = (0:max(unique_indices)-1)'; % 生成排名
sorted_ranks(end) = -1; % 最大值设为 -1

% 生成向量 l
l = sorted_ranks(unique_indices);
l = l';
% 找到 l 中的最大正整数 r
r = max(l(l > 0), [], 'omitnan'); % 如果没有正整数，返回空数组


disp('状态划分向量 :');
disp(l);
disp(['排序状态集从 0 到 r = ', num2str(r)]);

%% 构造alpha
t = l(i);
disp(['状态i所在的值集合序号为 ',num2str(t)])
alpha = R * double(l(i) ~= -1) + min(R,-2*x(i)) * double(t == -1);
disp(['这里的alpha = ',num2str(alpha)])


%% 得到状态集上的速率和x_com
% 提取 l 中的唯一值
all_unique = unique(l);

% 分离非负值并排序，最后添加 -1（如果存在）
non_negative = all_unique(all_unique ~= -1);
non_negative = sort(non_negative);
if any(all_unique == -1)
    new_unique = [non_negative, -1];
else
    new_unique = non_negative;
end

% 生成映射索引
[~, idx] = ismember(l, new_unique);

% 计算分组和
x_com = accumarray(idx', x');
x_com = x_com';
disp('原速率聚合向量:')
disp(x_com)


%% 最大值集合的最坏分布和
y = zeros(1,r+2);
y(end) = x_com(end) + alpha/2;
disp(['最大值集合的最坏分布和为 ',num2str(y(end))])

%% 最小值集合的最坏分布和
y(1) = posi(x_com(1) - alpha/2)*double(t ~= 0) + (x_com(1) - alpha/2)*double(t == 0);
disp(['最小值集合的最坏分布和为 ',num2str(y(1))])

%% 中间第s(s=1,2,...,r)个状态集合的最坏分布和
for s = 1:r
    if t == -1
        t = 10000;%如果i在最大状态集，那么本来我们规定他是-1，但为了下面的符号便捷，令其为一极大的数(这里取10000);
    end
    y(s+1) = x_com(s+1)*double(t < s)+...
        +(x_com(s+1)-posi(alpha/2 - sum(x_com(1:s))))*double(t == s)+...
        +posi(x_com(s+1)-posi(alpha/2 - sum(x_com(1:s))))*double(t > s);
    disp(['第 ',num2str(s),' 值集合的最坏分布和为 ',num2str(y(s+1))])
end

%% 输出最坏的分布(状态聚合形式)y和最大泛函值D
disp('最坏的概率分布y为(从最小状态集到最大状态集): ')
disp(y)
% 删除重复元素并排序
z = unique(k); % 删除重复元素并返回排序后的结果
disp('向量k聚合排序后为')
disp(z)
D = z*y';
disp(['最大的泛函值为: ',num2str(D)]);

end

function y = posi(x)
y = max(x,0);
end


