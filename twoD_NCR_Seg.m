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

function [traindata_para_cell, testdata_para_cell, Measure_Table, error_rate] = ...
    twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, SegPoint, Name)
%% 训练集和测试集初始化
[m1, ~] = size(traindata);             %训练集有m1组时间序列，每个时间序列有n维度；
[m2, n] = size(testdata);              %测试集有m2组时间序列，每个时间序列有n维度；

% disp(['输入数据集为：',Name]);
% disp(['训练集个数：',num2str(m1),'  测试集个数：',num2str(m2),'  数据维度为：',num2str(n)]);

traindata_diff = [zeros(m1, 1) diff(traindata, 1, 2)];     %差分训练集——一阶差分后的数据，最前端补零
testdata_diff = [zeros(m2, 1) diff(testdata, 1, 2)];       %差分测试集——一阶差分后的数据，最前端补零

%% 实现所有时序段的云模型特征变化

w = length(SegPoint) + 1; % w为分段数
SegPoint = [1 SegPoint n]; % 分段点坐标把头和尾补上

traindata_para_cell = cell(m1, 1);    %训练集云模型参数——每个元组内含有一个w*6的矩阵，代表一个时序数据的w段的云模型参数
testdata_para_cell  = cell(m2, 1);    %训练集云模型参数——每个元组内含有一个w*6的矩阵，代表一个时序数据的w段的云模型参数

for j = 1: m1 % 求训练集的第j个序列云参数
    Cloudparameter = zeros(w, 6);          %第j个序列的云参数，每一行是一段
    for k = 1: w % 求第k段的云参数
        x1 = traindata(j, SegPoint(k) + (k ~= 1): SegPoint(k + 1));       %x1代表训练集时序段的中间变量，第j个时间序列的第k段
        x2 = traindata_diff(j, SegPoint(k) + (k ~= 1): SegPoint(k + 1));  %x2代表差分训练集时序段的中间变量
        [Ex1, En1, He1, Ex2, En2, He2] = backward_twoD_Clouds(x1, x2);   %调用二维逆向云模型子函数
        Cloudparameter(k, :) = [Ex1, En1, He1, Ex2, En2, He2];
    end
    traindata_para_cell{j} = Cloudparameter;
end

for j = 1: m2 % 求测试集的第j个序列云参数
    Cloudparameter = zeros(w, 6);          %中间变量
    for k = 1: w
        x1 = testdata(j, SegPoint(k) + (k ~= 1): SegPoint(k + 1));       %x1代表训练集时序段的中间变量
        x2 =  testdata_diff(j, SegPoint(k) + (k ~= 1): SegPoint(k + 1));  %x2代表差分训练集时序段的中间变量
        [Ex1, En1, He1, Ex2, En2, He2] = backward_twoD_Clouds(x1, x2);   %调用二维逆向云模型子函数
        Cloudparameter(k,:) = [Ex1, En1, He1, Ex2, En2, He2];
    end
    testdata_para_cell{j} = Cloudparameter;
end

%% 实现两条时序数据云距离计算 
Measure_Table = zeros(m2, m1); %距离矩阵，第M行，N列，代表着第M个测试集时序数据与第N个训练集时序数据的重叠面积，值越大代表两组时序数据相似度越大。

for i=1: m2       %i代表测试集
    for j = 1: m1   %j代表训练集
        S = 0;     %S代表中间变量
        for k = 1: w
            C1_parameter = traindata_para_cell{j}(k,:);
            C2_parameter = testdata_para_cell{i}(k,:);
            [TOM] = Measure_for_twoD_Clouds(C1_parameter, C2_parameter);
            S = S + TOM;
        end
        Measure_Table(i, j) = sqrt(S / w) ;
    end
end

%% 实现1NN算法 
% K = 1; %1NN
[~, I] = sort(Measure_Table,2,'descend'); %每行相似度由大到小排序
Predictlabel = zeros(m2, 1);              %1NN算法的预测标签
for j = 1: m2
    Predictlabel(j, 1) = trainlabel(I(j, 1)) ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Step4  END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%%%%%%%%%%%%%%%%%%%%%%%% Step5：error rate计算 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    
    
    
    

