%{
���룺ѵ����--traindata ��ѵ������ǩ--trainlabel �����Լ�--testdata �����Լ���ǩ--testlabel�� �ֶ���--w��
        ���ݼ�����--Name���ַ��������ֶε�����--SegPoint(������)
������ֶκ��ѵ����Ԫ��--Cloudtraincell �� �ֶκ�Ĳ��Լ�Ԫ��--Cloudtestcell ��
        ���ƶȾ���--Measure_Table �� �����--error_rate��
ע��  ѵ����Ԫ��Ϊ m1*1����һάʱ�����ݵȷֳɶ�Σ�����һ������ÿһ�д���һ�Σ���w�У�����һ��Ԫ��cell�
        ���Լ�Ԫ��ͬѵ����Ԫ�飻
        ���ƶȾ���m2*m1�ľ��󣻵�i�У�j�У������ŵ�i�����Լ�ʱ���������j��ѵ����ʱ�����ݵ��ص������
        ֵԽ���������ʱ���������ƶ�Խ��
���ú�����backward_twoD_Clouds.m; Measure_for_twoD_Clouds.m;
%}

function [Measure_Table, error_rate] = ...
    twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, SegPoint, ~)
%% ѵ�����Ͳ��Լ���ʼ��
[m1, ~] = size(traindata);             % ѵ������m1��ʱ�����У�ÿ��ʱ��������nά�ȣ�
[m2, n] = size(testdata);              % ���Լ���m2��ʱ�����У�ÿ��ʱ��������nά�ȣ�

% disp(['�������ݼ�Ϊ��',Name]);
% disp(['ѵ����������',num2str(m1),'  ���Լ�������',num2str(m2),'  ����ά��Ϊ��',num2str(n)]);

traindata_diff = [zeros(m1, 1) diff(traindata, 1, 2)];     % ���ѵ��������һ�ײ�ֺ�����ݣ���ǰ�˲���
testdata_diff = [zeros(m2, 1) diff(testdata, 1, 2)];       % ��ֲ��Լ�����һ�ײ�ֺ�����ݣ���ǰ�˲���

%% ʵ������ʱ��ε���ģ�������仯��ѵ�����Ͳ��Լ���������������

w = length(SegPoint) + 1; % wΪ�ֶ���
SegPoint = [1 SegPoint n]; % �ֶε������ͷ��β����

train_Ex1 = zeros(m1, w); % ��ģ�Ͳ�����ʼ��
train_En1 = zeros(m1, w);
train_He1 = zeros(m1, w);
train_Ex2 = zeros(m1, w);
train_En2 = zeros(m1, w);
train_He2 = zeros(m1, w);

test_Ex1 = zeros(m2, w); % ��ģ�Ͳ�����ʼ��
test_En1 = zeros(m2, w);
test_He1 = zeros(m2, w);
test_Ex2 = zeros(m2, w);
test_En2 = zeros(m2, w);
test_He2 = zeros(m2, w);

for k = 1: w % ���k�ε��Ʋ���[Ex1, En1, He1, Ex2, En2, He2]���ֱ����traindata_para, testdata_para
    x1 = traindata(:, SegPoint(k) + (k ~= 1): SegPoint(k + 1));       % x1����ѵ����ʱ��ε��м��������j��ʱ�����еĵ�k��
    x2 = traindata_diff(:, SegPoint(k) + (k ~= 1): SegPoint(k + 1));  % x2������ѵ����ʱ��ε��м����
    [train_Ex1(:, k), train_En1(:, k), train_He1(:, k), train_Ex2(:, k), train_En2(:, k), train_He2(:, k)] = ...
        backward_twoD_Clouds(x1, x2);   % ���ö�ά������ģ���Ӻ���
    
    x1 = testdata(:, SegPoint(k) + (k ~= 1): SegPoint(k + 1));       % x1������Լ�ʱ��ε��м����
    x2 = testdata_diff(:, SegPoint(k) + (k ~= 1): SegPoint(k + 1));  % x2�����ֲ��Լ�ʱ��ε��м����
    [test_Ex1(:, k), test_En1(:, k), test_He1(:, k), test_Ex2(:, k), test_En2(:, k), test_He2(:, k)] = ...
        backward_twoD_Clouds(x1, x2);   % ���ö�ά������ģ���Ӻ���

end

%% ʵ������ʱ�������ƾ������ 
Measure_Table = zeros(m2, m1); % ������󣬵�M�У�N�У������ŵ�M�����Լ�ʱ���������N��ѵ����ʱ�����ݵ��ص������ֵԽ���������ʱ���������ƶ�Խ��

for i = 1: m2
    TOM = Measure_for_twoD_Clouds( ...
        train_Ex1, train_En1, train_He1, train_Ex2, train_En2, train_He2, ...
        repmat(test_Ex1(i, :), m1, 1), repmat(test_En1(i, :), m1, 1), repmat(test_He1(i, :), m1, 1), ...
            repmat(test_Ex2(i, :), m1, 1), repmat(test_En2(i, :), m1, 1), repmat(test_He2(i, :), m1, 1));
    Measure_Table(i, :) = sqrt(sum(TOM, 2)' / w);
end

%{
% for i = 1: m2       % i������Լ�
%     for j = 1: m1   % j����ѵ����
%         S = 0;     % S�����м����
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

%% ʵ��1NN�㷨 
% K = 1; %1NN
[~, I] = sort(Measure_Table, 2, 'descend'); % ÿ�����ƶ��ɴ�С����
Predictlabel = zeros(m2, 1);              % 1NN�㷨��Ԥ���ǩ
for j = 1: m2
    Predictlabel(j, 1) = trainlabel(I(j, 1)) ;
end

%% error rate����
correct_number = 0; %��ȷ����
for j = 1: m2
    if Predictlabel(j, 1) == testlabel(j, 1)
        correct_number = correct_number + 1 ;
    end
end
error_rate = (m2 - correct_number) / m2 ;

%% disp ��ʾ

% disp(['�ֶ��� w = ',num2str(w)]);
% disp(['correctNumber = ',num2str(correct_number)]);
% disp(['TestNumber = ',num2str(m2)]);
% disp(['correct rate = ',num2str(1-error_rate)]);
disp(['�ֶ��� ', num2str(w), ' error rate = ',num2str(error_rate)]);
    
end
    
    
    
    
    

