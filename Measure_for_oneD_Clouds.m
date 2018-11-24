%{
函数功能：计算一维云参数，可以并行计算
输入：六个参数矩阵，每一行是一个序列，许多行并行计算加快速度
输出：矩阵OD, u
%}

function [OD, u] = Measure_for_oneD_Clouds(Ex1, En1, He1, Ex2, En2, He2)
[m, w] = size(Ex1); % 获取并行计算的行数m，序列分段数w

En1(En1 == 0) = 10 ^-30;
En2(En2 == 0) = 10 ^-30;

SupC1 = Ex1 + 3 * En1;
InfC1 = Ex1 - 3 * En1;
SupC2 = Ex2 + 3 * En2;
InfC2 = Ex2 - 3 * En2;

OD = 2 * (min(SupC1, SupC2) - max(InfC1, InfC2)) ./ (6 * En1 + 6 * En2);

%% 
% 初始化x1, x2, n, u
x1 = zeros(m, w);
x2 = zeros(m, w);
n = zeros(m, w);
u = zeros(m, w);

a = zeros(m, w);
b = zeros(m, w);

index = En1 == En2; % 中间变量
    x1(index) = (Ex1(index) .* En2(index) + Ex2(index) .* En1(index)) ./ (En1(index) + En2(index));
    n(index) = 1; % 代表只有一个交点，当En1=En2时，才仅有一个根

index = En1 ~= En2; % 条件语句在并行运算里的改写
    a(index) = (Ex2(index) .* En1(index) - Ex1(index) .* En2(index)) ./ (En1(index) - En2(index));
    b(index) = (Ex1(index) .* En2(index) + Ex2(index) .* En1(index)) ./ (En1(index) + En2(index));
    x1(index) = min(a(index), b(index));
    x2(index) = max(a(index), b(index));   
    n(index) = 2; % 代表有2个交点

%%
% 只有一个交点时的u值
e = 10 ^-20; % 仅用来取消等号的判断
flag = (max(InfC1, InfC2) - e) < x1 & x1 < (min(SupC1, SupC2) + e); % 判断条件
index = n == 1 & flag; % 条件语句在并行运算里的改写
	u(index) = exp( -(x1(index) - Ex1(index)) .^2 ./ (2 * En1(index) .^2));
index = n == 1 & flag == 0;
	u(index) = 0;

% 有两个交点时的u值
u1 = zeros(m, w); % 初始化两个交点的u值
u2 = zeros(m, w);

index = n == 2 & flag; % 条件语句在并行运算里的改写
    u1(index) = exp( -(x1(index) - Ex1(index)) .^2 ./ (2 .* En1(index) .^2) );
flag = (max(InfC1, InfC2) - e) < x2 & x2 < (min(SupC1, SupC2) + e);
index = n == 2 & flag;
    u2(index) = exp( -(x2(index) - Ex2(index)) .^2 ./ (2 * En2(index) .^2) );
u1_u2 = max(u1, u2);
u(n == 2) = u1_u2(n == 2);

end




