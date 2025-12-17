function [ind,V,flag]=cal_index(x,F_func,opt)
% calculate index at x
g=F_func(x,opt);
ind=0;
V=[];
flag=1;

% 检查是否为临界点（容差可调整）
crit_tol = 1e-3;  % 保持原值
if norm(g) > crit_tol
    % disp('It is not a critical point')
    flag=0;
    return
end

%% 改进的参数设置
n = length(x);

% 方法1: 直接数值雅可比矩阵法（更可靠）
% 自适应有限差分步长
l = 1e-6 * max(1, norm(x));
if l == 0 || isnan(l)
    l = 1e-6;
end

% 计算雅可比矩阵
J = zeros(n, n);
for i = 1:n
    ei = zeros(n, 1);
    ei(i) = 1;
    % 中心差分
    F_plus = F_func(x + l * ei, opt);
    F_minus = F_func(x - l * ei, opt);
    J(:, i) = (F_plus - F_minus) / (2 * l);
end

% 计算特征值和特征向量
[V_eig, D] = eig(J);
lambda = diag(D);

% 计算指数（正实部特征值个数）
ind = sum(real(lambda) > 1e-8);

% 返回特征向量作为基
V = V_eig;

end