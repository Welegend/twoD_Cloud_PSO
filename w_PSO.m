% �������ܣ�������Ⱥ�㷨Ѱ��ʹ��λ��ģ���������С�ķֶε�
% ���룺ѵ����--traindata ��ѵ������ǩ--trainlabel �����Լ�--testdata �����Լ���ǩ--testlabel�� �ֶ���--w�� ���ݼ�����--Name���ַ�����
% ������������С�ķֶε�λ��������--SegPoint����С�������--error_rate
% ���ú�����twoD_NCR_Seg.m;

function [SegPoint, error_rate] = w_PSO(traindata, trainlabel, testdata, testlabel, w, Name)
dbstop if error % ����ʹ��
%% ����Ⱥ��ʼ��
[~, n] = size(traindata);

N = 20; % ��Ⱥ��ģ
D = w - 1; % ����ά�ȣ�w-1���㽫���зֳ�w��
T_init = 1;
T = 100;
Xmin = 1; % ��������ķ�Χ
Xmax = n - 1;
Vmin = -(n - 1);
Vmax = n - 1; % ���ӷ����ٶȵķ�Χ

% ww = 0.8; % ����ϵ��
c1 = 1.5; % ��֪ϵ��
c2 = 1.5; % ���ѧϰϵ��

% ��ʼ������Ⱥ��λ�ú��ٶ�
[~, popx] = sort(rand(N, n), 2);
popx = sort(popx(:, 1: D), 2); % ����Ⱥ�����ʼ����ÿһ�б�ʾÿ�����ӵķֶε����꣨�����ظ���
popv = round(rand(N, D) * (Vmax - Vmin) + Vmin); % ����Ⱥ�ٶȳ�ʼ����ÿһ�б�ʾÿ�����ӵ��ٶ�

% ÿ�����ӵ����λ�ó�ʼ��
pBest = popx;
pBestValue = zeros(N, 1); % ����ÿ�����ӵ���Ӧ�Ⱥ����������Ƕ�ά�Ƶ���
for i = 1: N
    [~, pBestValue(i)] = twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, popx(i, :), Name);
end

% ����Ⱥ��ȫ�����λ�ó�ʼ�����������������С�ķֶ����꣩
[gBestValue, index] = min(pBestValue);
gBest = pBest(index, :);

% ������������ŵ����ı仯���۲��Ƿ�������������
gBestValue_plot = zeros(T, 1);
figure;
set(gcf, 'unit', 'centimeters', 'Position', [10, 5, 9, 7]); % ����ͼ���С

%% ���¸����λ�ú��ٶ�
for t = 1: T % ��������T
    ww = 0.9 - (0.9 - 0.4) * T_init / T; % ����ϵ��,������������Ӷ��ݼ�
    
   %% �Ŵ��㷨��������λ�ú��ٶȣ���������Ⱥ�㷨֮ǰ��ʡȥ�˼���ظ��Ĳ���
    num = ceil(0.2 * N); % ��N������Ⱥ�г���ٷֱ�0.2�����ӣ���ĿΪnum������ȡ��
    %%%%%%%%%%%%%ѡ��%%%%%%%%%%%%%%%
    [~, I] = sort(pBestValue, 'descend'); % pBestValue��ı���Խ���ã�������num���ó�����popx_pool, popv_pool
    popx_pool = popx(I(1: num + 1), :); % num��һ����Ϊ���潻��׼����
%     popv_pool = popv(I(1: num + 1), :);

    %{
    % ���̶�ת��ʽ
    pBestValue_one = pBestValue / sum(pBestValue);
    pBestValue_one_cumsum = cumsum(pBestValue_one); % �����ۼƸ��ʣ�pBestValue��ĸ����ױ��鵽
    popx_pool = zeros(num, D);
    popv_pool = zeros(num, D);
    for j = 1: num % ���̶�תnum��
        k = find(pBestValue_one_cumsum >= rand, 1); % ��ȡ�±�ֵ
        popx_pool(j, :) = popx(k, :);
        popv_pool(j, :) = popv(k, :);
    end
    %}

    %%%%%%%%%%%%%����%%%%%%%%%%%%%%%
    for i = 1: num % ��popx_pool�ĵ�i�����ν���
        childx = [popx_pool(i, :), popx_pool(i + 1, :)];
        childx = childx(randperm(length(childx))); % ��popx_pool���д���
        childx = sort(reshape(childx, 2, []), 2); % �õ����н����ĺ��
        index = any(childx(:, 1: end - 1) == childx(:, 2: end), 2); % index��chidlx�����ظ�Ԫ�ص����е�����
        childx(index, :) = []; % ��ʱ��chidlxÿһ��Ԫ�ض��ǲ��ظ���
        x = size(childx, 1); % childx��ʣx��

        childx_error_rate = zeros(x, 1); % ���chidx��Ӧ�������
        for j = 1: x
            [~, childx_error_rate(j)] = twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, childx(j, :), Name);
        end
        popx_faAchi = [popx_pool(i, :); childx];
        error_rate_faAchi = [pBestValue(I(i)); childx_error_rate];

        % ���̶�ת��ʽ
        error_rate_pool = 1 ./ error_rate_faAchi; % ���̶�ת����Խ��Խ���׳鵽�������Ǹ�������С��
        error_rate_pool = error_rate_pool ./ sum(error_rate_pool);
        error_rate_cumsum = cumsum(error_rate_pool); % �����ۼƸ��ʣ�pBestValue��ĸ����ױ��鵽
        k = find(error_rate_cumsum >= rand, 1); % ��תһ����ȡ�±�ֵ
        popx(I(i), :) = popx_faAchi(k, :); % �������ӽ�����λ��
        if error_rate_faAchi(k) < pBestValue(I(i)) % ���¸�������
            pBest(I(i), :) = popx_faAchi(k, :);
            pBestValue(I(i)) = error_rate_faAchi(k);
            if pBestValue(I(i)) < gBestValue % ����ȫ������
                gBest = pBest(I(i), :);
                gBestValue = pBestValue(I(i));
            end
        end
        %{
        childv_flag = 0;
        cpoint = randperm(D - 1, 1); % ������ɽ����

        childx1 = sort([popx_pool(j, 1: cpoint), popx_pool(j + 1, cpoint + 1: end)]); % �����ĺ��
        if all(childx1(1: end - 1) ~= childx1(2: end)) % ����Ƿ����ظ�Ԫ��
            [~, childx1_error_rate] = twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, childx1, Name); % ���㽻���������
            if childx1_error_rate < pBestValue(I(j)) % ������������ʱ�С����ô���������Ľ��
                popx(I(j), :) = childx1;
                pBest(I(j), :) = childx1;
                pBestValue(I(j)) = childx1_error_rate;
                childv_flag = 1; % ������Ҫ����popv
            end
        end
        childx2 = sort([popx_pool(j + 1, 1: cpoint), popx_pool(j, cpoint + 1: end)]); % �����ĺ��
        if all(childx2(1: end - 1) ~= childx2(2: end)) % ����Ƿ����ظ�Ԫ��
            [~, childx2_error_rate] = twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, childx2, Name); % ���㽻���������
            if childx2_error_rate < pBestValue(I(j)) % ������������ʱ�С����ô���������Ľ��
                popx(I(j), :) = childx2;
                pBest(I(j), :) = childx2;
                pBestValue(I(j)) = childx2_error_rate;
                childv_flag = 1;
            end
        end
        %}
        if k ~= 1 % popx�仯��
            childv = ww * popx_faAchi(1, :) - popx_faAchi(k, :);
            popv(I(i), :) = childv;
            %{
            childv = (popv_pool(i, :) + popv_pool(i + 1, :)) .* norm(popv_pool(i, :), 2) ./ ...
                norm(popv_pool(i, :) + popv_pool(i + 1, :), 2); % ���ӷ�ĸ��̫�����׳���NaN
            if any(isnan(childv)) == 0 % ��ֹchildv�õ�NaN����
                popv(I(i), :) = round(childv);
            end
            %}
        end
    end
    
    %% ����Ⱥ�㷨���¸����λ�ú��ٶ�    
    for i = 1: N % һ��N�����ӣ��������
        flag = 1;
        while flag % Ϊ��ʵ��ֱ����ѭ��
            popv(i, :) = round(ww * popv(i, :) + c1 * rand * (pBest(i, :) - popx(i, :)) + c2 * rand * (gBest - popx(i, :)));
            popv(i, popv(i, :) > Vmax) = randi(Vmax, 1, sum(popv(i, :) > Vmax)); % ����û����һ��������ѭ����
            popv(i, popv(i, :) < Vmin) = -randi(Vmax, 1, sum(popv(i, :) < Vmin));
            
            popx(i, :) = sort(popx(i, :) + popv(i, :));
            popx(i, popx(i, :) > Xmax) = Xmax;
            popx(i, popx(i, :) < Xmin) = Xmin;
            
            if w == 2 % ��ʱpopxΪ��������һ�����ظ�
                flag = 0;
            elseif all(popx(i, 1: end - 1) ~= popx(i, 2: end)) % popx��i����������Ԫ�ض����ظ�
                flag = 0;
            end
        end
        
        % ����������ʷ����
        disp(['���µ����� ', num2str(i)]);% ����ʱʹ��
        
        [~, error_rate] = twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, popx(i, :), Name);
        if error_rate < pBestValue(i) % ���¸�������
            pBest(i, :) = popx(i, :);
            pBestValue(i) = error_rate;
            if pBestValue(i) < gBestValue % ����ȫ�����ţ����ﲻ�ܰ�pBestValue(i)����error_rate����Ϊǰ�潻�����ı���pBestValue
                gBest = pBest(i, :);
                gBestValue = pBestValue(i);
            end
        end
    end
    T_init = T_init + 1;
    gBestValue_plot(t) = gBestValue;
    
    %% ����ʱʹ��
    disp(['�ֶ��� ', num2str(w), ' �������� = ',num2str(t), ' ȫ������ֵ = ',num2str(gBestValue)]);

    plot(gBestValue_plot);
    xlabel('��������');
    ylabel('�����');
    xlim([1 T]); % ����������
    title(['���ݼ� ', Name, ' �ֶ��� ', num2str(w)]); % ͼ��
    set(gca, 'FontName', '����', 'FontSize', 10.5);
    drawnow
end

%% ����Ⱥ������ɣ����ȫ������ֵ
SegPoint = gBest;
error_rate = gBestValue;

end