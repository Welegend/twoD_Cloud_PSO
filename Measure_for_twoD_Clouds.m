%{ 
调用函数：Measure_for_oneD_Clouds.m;
输入：6个训练集参数，6个测试集参数，都是矩阵，行数是数据集序列个数，列数是六个参数
输出：矩阵TOM
%}

function TOM = Measure_for_twoD_Clouds( ...
    train_Ex1, train_En1, train_He1, train_Ex2, train_En2, train_He2, ...
    test_Ex1, test_En1, test_He1, test_Ex2, test_En2, test_He2)

% 注意输入参数是训练集和测试集两朵交叉的
[OD1, u1] = Measure_for_oneD_Clouds(train_Ex1, train_En1, train_He1, test_Ex1, test_En1, test_He1);
[OD2, u2] = Measure_for_oneD_Clouds(train_Ex2, train_En2, train_He2, test_Ex2, test_En2, test_He2);

% 输出参数
OM1 = OD1 .* u1;
OM2 = OD2 .* u2;
TOM = OM1 .* OM2;

end