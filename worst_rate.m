function [y,D,n,i,r] = worst_rate(x,k,R)

%% Ԥ����
n = length(x);
disp(['��� ',num2str(n),' ά�������Ϊ'])
disp(x)
i = find(x < 0);
disp(['��ǰ״̬����iΪ ',num2str(i)]);
disp(['TVԼ������R(i)Ϊ ',num2str(R)]);
disp(['��� ',num2str(n),' άϵ������Ϊ'])
disp(k)

%% ����״̬�� 
% ȷ�������ǷǸ�����
if any(k < 0)
    error('�������� k �����ǷǸ���');
end

% ��ȡΨһֵ������
[~, ~, unique_indices] = unique(k); % ��ȡΨһֵ����
sorted_ranks = (0:max(unique_indices)-1)'; % ��������
sorted_ranks(end) = -1; % ���ֵ��Ϊ -1

% �������� l
l = sorted_ranks(unique_indices);
l = l';
% �ҵ� l �е���������� r
r = max(l(l > 0), [], 'omitnan'); % ���û�������������ؿ�����


disp('״̬�������� :');
disp(l);
disp(['����״̬���� 0 �� r = ', num2str(r)]);

%% ����alpha
t = l(i);
disp(['״̬i���ڵ�ֵ�������Ϊ ',num2str(t)])
alpha = R * double(l(i) ~= -1) + min(R,-2*x(i)) * double(t == -1);
disp(['�����alpha = ',num2str(alpha)])


%% �õ�״̬���ϵ����ʺ�x_com
% ��ȡ l �е�Ψһֵ
all_unique = unique(l);

% ����Ǹ�ֵ������������ -1��������ڣ�
non_negative = all_unique(all_unique ~= -1);
non_negative = sort(non_negative);
if any(all_unique == -1)
    new_unique = [non_negative, -1];
else
    new_unique = non_negative;
end

% ����ӳ������
[~, idx] = ismember(l, new_unique);

% ��������
x_com = accumarray(idx', x');
x_com = x_com';
disp('ԭ���ʾۺ�����:')
disp(x_com)


%% ���ֵ���ϵ���ֲ���
y = zeros(1,r+2);
y(end) = x_com(end) + alpha/2;
disp(['���ֵ���ϵ���ֲ���Ϊ ',num2str(y(end))])

%% ��Сֵ���ϵ���ֲ���
y(1) = posi(x_com(1) - alpha/2)*double(t ~= 0) + (x_com(1) - alpha/2)*double(t == 0);
disp(['��Сֵ���ϵ���ֲ���Ϊ ',num2str(y(1))])

%% �м��s(s=1,2,...,r)��״̬���ϵ���ֲ���
for s = 1:r
    if t == -1
        t = 10000;%���i�����״̬������ô�������ǹ涨����-1����Ϊ������ķ��ű�ݣ�����Ϊһ�������(����ȡ10000);
    end
    y(s+1) = x_com(s+1)*double(t < s)+...
        +(x_com(s+1)-posi(alpha/2 - sum(x_com(1:s))))*double(t == s)+...
        +posi(x_com(s+1)-posi(alpha/2 - sum(x_com(1:s))))*double(t > s);
    disp(['�� ',num2str(s),' ֵ���ϵ���ֲ���Ϊ ',num2str(y(s+1))])
end

%% �����ķֲ�(״̬�ۺ���ʽ)y����󷺺�ֵD
disp('��ĸ��ʷֲ�yΪ(����С״̬�������״̬��): ')
disp(y)
% ɾ���ظ�Ԫ�ز�����
z = unique(k); % ɾ���ظ�Ԫ�ز����������Ľ��
disp('����k�ۺ������Ϊ')
disp(z)
D = z*y';
disp(['���ķ���ֵΪ: ',num2str(D)]);

end

function y = posi(x)
y = max(x,0);
end


