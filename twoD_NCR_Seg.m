%{
输入：训练集--traindata ；训练集标签--trainlabel ；测试集--testdata ；测试集标签--testlabel； 分段数--w；
        数据集名称--Name（字符串）；分段点坐标--SegPoint(行向量)
输出：分段后的训练集元组--Cloudtraincell ； 分段后的测试集元组--Cloudtestcell ；
        相似度矩阵--Measure_Table ； 误差率--error_rate；
注：  训练集元组为 m1*1，将一维时序数据等分成多段，构成一个矩阵（每一行代表一段，共w行）存在一个元组cell里。
        测试集元组同训练集元组；
        相似度矩阵：m2*m1的矩阵；第i行，j列，代表着第i个测试集时序数据与第j个训练集时序数据的重叠面积，
        值越大代表两组时序数据相似度越大
调用函数：backward_twoD_Clouds.m; Measure_for_twoD_Clouds.m;
%}

function [Measure_Table, error_rate] = ...
    twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, SegPoint, ~)
%% 训练集和测试集初始化
[m1, ~] = size(traindata);             % 训练集有m1组时间序列，每个时间序列有n维度；
[m2, n] = size(testdata);              % 测试集有m2组时间序列，每个时间序列有n维度；

% disp(['输入数据集为：',Name]);
% disp(['训练集个数：',num2str(m1),'  测试集个数：',num2str(m2),'  数据维度为：',num2str(n)]);

traindata_diff = [zeros(m1, 1) diff(traindata, 1, 2)];     % 差分训练集——一阶差分后的数据，最前端补零
testdata_diff = [zeros(m2, 1) diff(testdata, 1, 2)];       % 差分测试集——一阶差分后的数据，最前端补零

%% 实现所有时序段的云模型特征变化，训练集和测试集的所有六个参数

w = length(SegPoint) + 1; % w为分段数
SegPoint = [1 SegPoint n]; % 分段点坐标把头和尾补上

train_Ex1 = zeros(m1, w); % 云模型参数初始化
train_En1 = zeros(m1, w);
train_He1 = zeros(m1, w);
train_Ex2 = zeros(m1, w);
train_En2 = zeros(m1, w);
train_He2 = zeros(m1, w);

test_Ex1 = zeros(m2, w); % 云模型参数初始化
test_En1 = zeros(m2, w);
test_He1 = zeros(m2, w);
test_Ex2 = zeros(m2, w);
test_En2 = zeros(m2, w);
test_He2 = zeros(m2, w);

for k = 1: w % 求第k段的云参数[Ex1, En1, He1, Ex2, En2, He2]，分别放在traindata_para, testdata_para
    x1 = traindata(:, SegPoint(k) + (k ~= 1): SegPoint(k + 1));       % x1代表训练集时序段的中间变量，第j个时间序列的第k段
    x2 = traindata_diff(:, SegPoint(k) + (k ~= 1): SegPoint(k + 1));  % x2代表差分训练集时序段的中间变量
    [train_Ex1(:, k), train_En1(:, k), train_He1(:, k), train_Ex2(:, k), train_En2(:, k), train_He2(:, k)] = ...
        backward_twoD_Clouds(x1, x2);   % 调用二维逆向云模型子函数
    
    x1 = testdata(:, SegPoint(k) + (k ~= 1): SegPoint(k + 1));       % x1代表测试集时序段的中间变量
    x2 = testdata_diff(:, SegPoint(k) + (k ~= 1): SegPoint(k + 1));  % x2代表差分测试集时序段的中间变量
    [test_Ex1(:, k), test_En1(:, k), test_He1(:, k), test_Ex2(:, k), test_En2(:, k), test_He2(:, k)] = ...
        backward_twoD_Clouds(x1, x2);   % 调用二维逆向云模型子函数

end

%% 实现两条时序数据云距离计算 
Measure_Table = zeros(m2, m1); % 距离矩阵，第M行，N列，代表着第M个测试集时序数据与第N个训练集时序数据的重叠面积，值越大代表两组时序数据相似度越大。

for i = 1: m2
    TOM = Measure_for_twoD_Clouds( ...
        train_Ex1, train_En1, train_He1, train_Ex2, train_En2, train_He2, ...
        repmat(test_Ex1(i, :), m1, 1), repmat(test_En1(i, :), m1, 1), repmat(test_He1(i, :), m1, 1), ...
            repmat(test_Ex2(i, :), m1, 1), repmat(test_En2(i, :), m1, 1), repmat(test_He2(i, :), m1, 1));
    Measure_Table(i, :) = sqrt(sum(TOM, 2)' / w);
end

%{
% for i = 1: m2       % i代表测试集
%     for j = 1: m1   % j代表训练集
%         S = 0;     % S代表中间变量
%         for k = 1: w
%             C1_parameter = traindata_para(j,  k * 6 - 5: k * 6);
%             C2_parameter = testdata_para(i,  k * 6 - 5: k * 6);
%             [TOM] = Measure_for_twoD_Clouds(C1_parameter, C2_parameter);
%             S = S + TOM;
%         end
%         Measure_Table(i, j) = sqrt(S / w) ;
%     end
% end
%}

%% 实现1NN算法 
% K = 1; %1NN
[~, I] = sort(Measure_Table, 2, 'descend'); % 每行相似度由大到小排序
Predictlabel = zeros(m2, 1);              % 1NN算法的预测标签
for j = 1: m2
    Predictlabel(j, 1) = trainlabel(I(j, 1)) ;
end

%% error rate计算
correct_number = 0; %正确个数
for j = 1: m2
    if Predictlabel(j, 1) == testlabel(j, 1)
        correct_number = correct_number + 1 ;
    end
end
error_rate = (m2 - correct_number) / m2 ;

%% disp 显示

% disp(['分段数 w = ',num2str(w)]);
% disp(['correctNumber = ',num2str(correct_number)]);
% disp(['TestNumber = ',num2str(m2)]);
% disp(['correct rate = ',num2str(1-error_rate)]);
disp(['分段数 ', num2str(w), ' error rate = ',num2str(error_rate)]);
    
end
    
    
    
    
    

