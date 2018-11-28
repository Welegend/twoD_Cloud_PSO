% 函数功能：用粒子群算法寻找使二位云模型误差率最小的分段点
% 输入：训练集--traindata ；训练集标签--trainlabel ；测试集--testdata ；测试集标签--testlabel； 分段数--w； 数据集名称--Name（字符串）
% 输出：误差率最小的分段点位置行向量--SegPoint；最小的误差率--error_rate
% 调用函数：twoD_NCR_Seg.m;

function [SegPoint, error_rate] = w_PSO(traindata, trainlabel, testdata, testlabel, w, Name)
dbstop if error % 调试使用
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
    [~, pBestValue(i)] = twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, popx(i, :), Name);
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
    
   %% 遗传算法更新粒子位置和速度，放在粒子群算法之前，省去了检测重复的步骤
    r1 = rand();
    num = ceil(0.2 * N); % 从N个粒子群中抽出百分比0.2的粒子，数目为num，向上取整
    if r1 < 0.6 % 变异概率
        %%%%%%%%%%%%%选择%%%%%%%%%%%%%%%
        [~, I] = sort(pBestValue, 'descend'); % pBestValue大的表现越不好，把最大的num个拿出来给popx_pool, popv_pool
        popx_pool = popx(I(1: num + 1), :); % num多一个是为下面交叉准备的
        popv_pool = popv(I(1: num + 1), :);

        %{
        % 轮盘赌转方式
        pBestValue_one = pBestValue / sum(pBestValue);
        pBestValue_one_cumsum = cumsum(pBestValue_one); % 计算累计概率，pBestValue大的更容易被抽到
        popx_pool = zeros(num, D);
        popv_pool = zeros(num, D);
        for j = 1: num % 轮盘赌转num次
            k = find(pBestValue_one_cumsum >= rand, 1); % 提取下标值
            popx_pool(j, :) = popx(k, :);
            popv_pool(j, :) = popv(k, :);
        end
        %}

        %%%%%%%%%%%%%交叉变异%%%%%%%%%%%%%%%
        for j = 1: num
            childv_flag = 0;
            
            childx = [popx_pool(j, :), popx_pool(j + 1, :)];
            childx = childx(randperm(length(childx))); % 把popx_pool两行打乱
            childx1 = sort(childx(1: D)); % 交叉后的后代
            if all(childx1(1: end - 1) ~= childx1(2: end)) % 检测是否有重复元素
                [~, childx1_error_rate] = twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, childx1, Name); % 计算交叉后的误差率
                if childx1_error_rate < pBestValue(I(j)) % 如果交叉后误差率变小，那么保留交叉后的结果
                    popx(I(j), :) = childx1;
                    pBest(I(j), :) = childx1;
                    pBestValue(I(j)) = childx1_error_rate;
                    childv_flag = 1; % 后面需要更新popv
                end
            end
            childx2 = sort(childx(D + 1: end)); % 交叉后的后代
            if all(childx2(1: end - 1) ~= childx2(2: end)) % 检测是否有重复元素
                [~, childx2_error_rate] = twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, childx2, Name); % 计算交叉后的误差率
                if childx2_error_rate < pBestValue(I(j)) % 如果交叉后误差率变小，那么保留交叉后的结果
                    popx(I(j), :) = childx2;
                    pBest(I(j), :) = childx2;
                    pBestValue(I(j)) = childx2_error_rate;
                    childv_flag = 1;
                end
            end
            %{
            childv_flag = 0;
            cpoint = randperm(D - 1, 1); % 随机生成交叉点
            
            childx1 = sort([popx_pool(j, 1: cpoint), popx_pool(j + 1, cpoint + 1: end)]); % 交叉后的后代
            if all(childx1(1: end - 1) ~= childx1(2: end)) % 检测是否有重复元素
                [~, childx1_error_rate] = twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, childx1, Name); % 计算交叉后的误差率
                if childx1_error_rate < pBestValue(I(j)) % 如果交叉后误差率变小，那么保留交叉后的结果
                    popx(I(j), :) = childx1;
                    pBest(I(j), :) = childx1;
                    pBestValue(I(j)) = childx1_error_rate;
                    childv_flag = 1; % 后面需要更新popv
                end
            end
            childx2 = sort([popx_pool(j + 1, 1: cpoint), popx_pool(j, cpoint + 1: end)]); % 交叉后的后代
            if all(childx2(1: end - 1) ~= childx2(2: end)) % 检测是否有重复元素
                [~, childx2_error_rate] = twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, childx2, Name); % 计算交叉后的误差率
                if childx2_error_rate < pBestValue(I(j)) % 如果交叉后误差率变小，那么保留交叉后的结果
                    popx(I(j), :) = childx2;
                    pBest(I(j), :) = childx2;
                    pBestValue(I(j)) = childx2_error_rate;
                    childv_flag = 1;
                end
            end
            %}
            if childv_flag
                childv = (popv_pool(j, :) + popv_pool(j + 1, :)) .* norm(popv_pool(j, :), 2) ./ ...
                    norm(popv_pool(j, :) + popv_pool(j + 1, :), 2); % 分子分母数太大容易出现NaN
                if any(isnan(childv)) == 0 % 防止childv得到NaN出错
                    popv(I(j), :) = round(childv);
                end
            end
        end
    end
    
    %% 粒子群算法更新个体的位置和速度    
    for i = 1: N % 一共N个粒子，逐个更新
        flag = 1;
        while flag % 为了实现直到型循环
            popv(i, :) = round(ww * popv(i, :) + c1 * rand * (pBest(i, :) - popx(i, :)) + c2 * rand * (gBest - popx(i, :)));
            popv(i, popv(i, :) > Vmax) = randi(Vmax, 1, sum(popv(i, :) > Vmax)); % 这里没加有一次跳不出循环了
            popv(i, popv(i, :) < Vmin) = -randi(Vmax, 1, sum(popv(i, :) < Vmin));
            
            popx(i, :) = sort(popx(i, :) + popv(i, :));
            popx(i, popx(i, :) > Xmax) = Xmax;
            popx(i, popx(i, :) < Xmin) = Xmin;
            
            if w == 2 % 此时popx为列向量，一定不重复
                flag = 0;
            elseif all(popx(i, 1: end - 1) ~= popx(i, 2: end)) % popx第i行所有相邻元素都不重复
                flag = 0;
            end
        end
        
        % 更新粒子历史最优
        disp(['更新到粒子 ', num2str(i)]);% 调试时使用
        
        [~, error_rate] = twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, popx(i, :), Name);
        if error_rate < pBestValue(i) % 更新个体最优
            pBest(i, :) = popx(i, :);
            pBestValue(i) = error_rate;
            if pBestValue(i) < gBestValue % 更新全局最优，这里不能把pBestValue(i)换成error_rate，因为前面交叉变异改变了pBestValue
                gBest = pBest(i, :);
                gBestValue = pBestValue(i);
            end
        end
    end
    T_init = T_init + 1;
    gBestValue_plot(t) = gBestValue;
    
    %% 调试时使用
    disp(['分段数 ', num2str(w), ' 迭代次数 = ',num2str(t), ' 全局最优值 = ',num2str(gBestValue)]);

    plot(gBestValue_plot);
    xlabel('迭代次数');
    ylabel('误差率');
    xlim([1 T]); % 坐标轴设置
    title(['数据集 ', Name, ' 分段数 ', num2str(w)]); % 图题
    set(gca, 'FontName', '宋体', 'FontSize', 10.5);
    drawnow
end

%% 粒子群迭代完成，输出全局最优值
SegPoint = gBest;
error_rate = gBestValue;

end