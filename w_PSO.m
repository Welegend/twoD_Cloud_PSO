% �������ܣ�������Ⱥ�㷨Ѱ��ʹ��λ��ģ���������С�ķֶε�
% ���룺ѵ����--traindata ��ѵ������ǩ--trainlabel �����Լ�--testdata �����Լ���ǩ--testlabel�� �ֶ���--w�� ���ݼ�����--Name���ַ�����
% ������������С�ķֶε�λ��������--SegPoint����С�������--error_rate
% ���ú�����twoD_NCR_PSO.m;

function [SegPoint, error_rate] = w_PSO(traindata, trainlabel, testdata, testlabel, w, Name)
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
    [~, ~, ~, pBestValue(i)] = twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, popx(i, :), Name);
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
    for i = 1: N % һ��N�����ӣ��������
        
        % ���¸����λ�ú��ٶ�
        flag = 1;
        while flag % Ϊ��ʵ��ֱ����ѭ��
            popv(i, :) = round(ww * popv(i, :) + c1 * rand * (pBest(i, :) - popx(i, :)) + c2 * rand * (gBest - popx(i, :)));
            popv(i, popv(i, :) > Vmax) = randperm(Vmax, sum(popv(i, :) > Vmax)); % ����û����һ��������ѭ����
            popv(i, popv(i, :) < Vmin) = -randperm(Vmax, sum(popv(i, :) < Vmin));
            
            popx(i, :) = sort(popx(i, :) + popv(i, :));
            popx(i, popx(i, :) > Xmax) = Xmax;
            popx(i, popx(i, :) < Xmin) = Xmin;
            
            if w == 2 % ��ʱpopxΪ��������һ�����ظ�
                flag = 0;
            elseif all(popx(:, 1: end - 1) ~= popx(:, 2: end)) % ��������Ԫ�ض����ظ�����popx��������
                flag = 0;
            end
        end

        % ����������ʷ����
        disp(['���µ����� ', num2str(i)]);
        [~, ~, ~, error_rate] = twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, popx(i, :), Name);
        if error_rate < pBestValue(i) % ���¸�������
            pBest(i, :) = popx(i, :);
            pBestValue(i) = error_rate;
            if error_rate < gBestValue % ����ȫ������
                gBest = pBest(i, :);
                gBestValue = error_rate;
            end
        end
    end
    T_init = T_init + 1;
    gBestValue_plot(t) = gBestValue;
    
    % ����ʱʹ��
    disp(['�ֶ��� ', num2str(w), ' �������� = ',num2str(t), ' ȫ������ֵ = ',num2str(gBestValue)]);
    
    plot(gBestValue_plot);
    xlabel('��������');
    ylabel('�����');
    xlim([1 T]); % ����������
    title(['���ݼ� ', Name, ' �ֶ��� ', num2str(w)]); % ͼ��
    set(gca, 'FontName', '����', 'FontSize', 10.5);
end

%% ����Ⱥ������ɣ����ȫ������ֵ
SegPoint = gBest;
error_rate = gBestValue;

end