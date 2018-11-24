%{
�������ܣ�����һά�Ʋ��������Բ��м���
���룺������������ÿһ����һ�����У�����в��м���ӿ��ٶ�
���������OD, u
%}

function [OD, u] = Measure_for_oneD_Clouds(Ex1, En1, He1, Ex2, En2, He2)
[m, w] = size(Ex1); % ��ȡ���м��������m�����зֶ���w

En1(En1 == 0) = 10 ^-30;
En2(En2 == 0) = 10 ^-30;

SupC1 = Ex1 + 3 * En1;
InfC1 = Ex1 - 3 * En1;
SupC2 = Ex2 + 3 * En2;
InfC2 = Ex2 - 3 * En2;

OD = 2 * (min(SupC1, SupC2) - max(InfC1, InfC2)) ./ (6 * En1 + 6 * En2);

%% 
% ��ʼ��x1, x2, n, u
x1 = zeros(m, w);
x2 = zeros(m, w);
n = zeros(m, w);
u = zeros(m, w);

a = zeros(m, w);
b = zeros(m, w);

index = En1 == En2; % �м����
    x1(index) = (Ex1(index) .* En2(index) + Ex2(index) .* En1(index)) ./ (En1(index) + En2(index));
    n(index) = 1; % ����ֻ��һ�����㣬��En1=En2ʱ���Ž���һ����

index = En1 ~= En2; % ��������ڲ���������ĸ�д
    a(index) = (Ex2(index) .* En1(index) - Ex1(index) .* En2(index)) ./ (En1(index) - En2(index));
    b(index) = (Ex1(index) .* En2(index) + Ex2(index) .* En1(index)) ./ (En1(index) + En2(index));
    x1(index) = min(a(index), b(index));
    x2(index) = max(a(index), b(index));   
    n(index) = 2; % ������2������

%%
% ֻ��һ������ʱ��uֵ
e = 10 ^-20; % ������ȡ���Ⱥŵ��ж�
flag = (max(InfC1, InfC2) - e) < x1 & x1 < (min(SupC1, SupC2) + e); % �ж�����
index = n == 1 & flag; % ��������ڲ���������ĸ�д
	u(index) = exp( -(x1(index) - Ex1(index)) .^2 ./ (2 * En1(index) .^2));
index = n == 1 & flag == 0;
	u(index) = 0;

% ����������ʱ��uֵ
u1 = zeros(m, w); % ��ʼ�����������uֵ
u2 = zeros(m, w);

index = n == 2 & flag; % ��������ڲ���������ĸ�д
    u1(index) = exp( -(x1(index) - Ex1(index)) .^2 ./ (2 .* En1(index) .^2) );
flag = (max(InfC1, InfC2) - e) < x2 & x2 < (min(SupC1, SupC2) + e);
index = n == 2 & flag;
    u2(index) = exp( -(x2(index) - Ex2(index)) .^2 ./ (2 * En2(index) .^2) );
u1_u2 = max(u1, u2);
u(n == 2) = u1_u2(n == 2);

end




