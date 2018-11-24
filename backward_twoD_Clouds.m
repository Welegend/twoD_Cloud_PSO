%{
函数功能：输入原始数据集和差分数据集的每一段，二维逆向云模型计算参数，可以并行计算数据集所有的序列
输入：x1的第i行为原始数据集第i个序列的某一段，x2的第i行为差分数据集第i个序列的某一段
输出：列向量 Ex1, En1, He1, Ex2, En2, He2
%}
function [Ex1, En1, He1, Ex2, En2, He2] = backward_twoD_Clouds(x1, x2)

 % 输入两组等长时间序列段
 Ex1 = mean(x1, 2); % 按行（序列）取x1的平均值
 Ex2 = mean(x2, 2);
 En1 = mean(abs(x1 - Ex1), 2) * sqrt(pi / 2);
 En2 = mean(abs(x2 - Ex2), 2) * sqrt(pi / 2);
 He1 = sqrt(var(x1, 0, 2) - En1 .^2);
 He2 = sqrt(var(x2, 0, 2) - En2 .^2);
 
end