% 函数功能：用粒子群算法寻找使二位云模型误差率最小的分段点
% 输入：训练集--traindata ；训练集标签--trainlabel ；测试集--testdata ；测试集标签--testlabel； 分段数--w； 数据集名称--Name（字符串）
% 输出：误差率最小的分段点位置行向量--SegPoint；最小的误差率--error_rate
% 调用函数：twoD_NCR_PSO.m;

function [SegPoint, error_rate] = w_PSO(traindata, trainlabel, testdata, testlabel, w, Name)
%% 粒子群初始化
[~, n] = size(traindata);

N = 20; % 种群规模
D = w - 1; % 粒子维度，w-1个点将序列分成w段
T_init = 1;
T = 100;
Xmin = 1; % 粒子坐标的范围
Xmax = n - 1;
Vmin = -(n - 1);
Vmax = n - 1; % 粒子飞行速度的范围

% ww = 0.8; % 惯性系数
c1 = 1.5; % 认知系数
c2 = 1.5; % 社会学习系数

% 初始化粒子群的位置和速度
[~, popx] = sort(rand(N, n), 2);
popx = sort(popx(:, 1: D), 2); % 粒子群坐标初始化，每一行表示每个粒子的分段点坐标（不能重复）
popv = round(rand(N, D) * (Vmax - Vmin) + Vmin); % 粒子群速度初始化，每一行表示每个粒子的速度

% 每个粒子的最佳位置初始化
pBest = popx;
pBestValue = zeros(N, 1); % 计算每个粒子的适应度函数（这里是二维云的误差）
for i = 1: N
    [~, ~, ~, pBestValue(i)] = twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, popx(i, :), Name);
end

% 粒子群的全局最佳位置初始化（这里是误差率最小的分段坐标）
[gBestValue, index] = min(pBestValue);
gBest = pBest(index, :);

% 绘制误差率随着迭代的变化，观察是否收敛，调试用
gBestValue_plot = zeros(T, 1);
figure;
set(gcf, 'unit', 'centimeters', 'Position', [10, 5, 9, 7]); % 设置图像大小

%% 更新个体的位置和速度
for t = 1: T % 迭代次数T
    ww = 0.9 - (0.9 - 0.4) * T_init / T; % 惯性系数,随迭代次数增加而递减
    for i = 1: N % 一共N个粒子，逐个更新
        
        % 更新个体的位置和速度
        flag = 1;
        while flag % 为了实现直到型循环
            popv(i, :) = round(ww * popv(i, :) + c1 * rand * (pBest(i, :) - popx(i, :)) + c2 * rand * (gBest - popx(i, :)));
            popv(i, popv(i, :) > Vmax) = randperm(Vmax, sum(popv(i, :) > Vmax)); % 这里没加有一次跳不出循环了
            popv(i, popv(i, :) < Vmin) = -randperm(Vmax, sum(popv(i, :) < Vmin));
            
            popx(i, :) = sort(popx(i, :) + popv(i, :));
            popx(i, popx(i, :) > Xmax) = Xmax;
            popx(i, popx(i, :) < Xmin) = Xmin;
            
            if w == 2 % 此时popx为列向量，一定不重复
                flag = 0;
            elseif all(popx(:, 1: end - 1) ~= popx(:, 2: end)) % 所有相邻元素都不重复，或popx是列向量
                flag = 0;
            end
        end

        % 更新粒子历史最优
        disp(['更新到粒子 ', num2str(i)]);
        [~, ~, ~, error_rate] = twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, popx(i, :), Name);
        if error_rate < pBestValue(i) % 更新个体最优
            pBest(i, :) = popx(i, :);
            pBestValue(i) = error_rate;
            if error_rate < gBestValue % 更新全局最优
                gBest = pBest(i, :);
                gBestValue = error_rate;
            end
        end
    end
    T_init = T_init + 1;
    gBestValue_plot(t) = gBestValue;
    
    % 调试时使用
    disp(['分段数 ', num2str(w), ' 迭代次数 = ',num2str(t), ' 全局最优值 = ',num2str(gBestValue)]);
    
    plot(gBestValue_plot);
    xlabel('迭代次数');
    ylabel('误差率');
    xlim([1 T]); % 坐标轴设置
    title(['数据集 ', Name, ' 分段数 ', num2str(w)]); % 图题
    set(gca, 'FontName', '宋体', 'FontSize', 10.5);
end

%% 粒子群迭代完成，输出全局最优值
SegPoint = gBest;
error_rate = gBestValue;

end