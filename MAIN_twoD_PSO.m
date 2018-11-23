% 函数功能：用PSO算法计算分段数2-9时，最佳分段点，一次只能测试其中的一个数据集
% 输入：数据集名称元胞数组
% 输出：对应分段数w下，最佳分段点SegPoint_cell，最小误差率error_rate，并画出w下误差率迭代曲线

function [w, SegPoint_cell, error_rate] = MAIN_twoD_PSO(newName_Dataset)
%% 加载数据集
i = 1; % 加载newName_Dataset中的第一个数据集
Name = newName_Dataset{i, 2}; % 数据集名称
Path = 'F:\研究生\数据挖掘课题组\实验模型\二维云模型\UCR_TS_Archive_2015';
eval(['load ', Path, '\', newName_Dataset{i, 2}, '\',Name, '_TRAIN']);
eval(['load ', Path, '\', newName_Dataset{i, 2}, '\',Name, '_TEST']);

%% 由数据集得到w_PSO的输入变量
Rowtraindata = eval([newName_Dataset{i, 2}, '_TRAIN']);
Rowtestdata = eval([newName_Dataset{i, 2}, '_TEST']);

Rowtraindata = sortrows(Rowtraindata, 1) ;             %根据类排布在一起（按第一列对行进行重排列）
Rowtestdata = sortrows(Rowtestdata, 1) ;

traindata = Rowtraindata(:, 2: end);       %训练集——每行代表一组时序数据
testdata = Rowtestdata(:, 2: end);         %测试集——每行代表一组时序数据
trainlabel = Rowtraindata(:, 1);          %训练集标签
testlabel = Rowtestdata(:, 1);            %测试集标签

%% 用PSO算法求分段数由2到9的最佳分段点和误差率
w = ([4, 8])'; % w为要寻优的分段数
SegPoint_cell = cell(length(w), 1);
error_rate = zeros(length(w), 1);
for i = 1: length(w)
    [SegPoint_cell{i}, error_rate(i)] = w_PSO(traindata, trainlabel, testdata, testlabel, w(i), Name);
end

end




