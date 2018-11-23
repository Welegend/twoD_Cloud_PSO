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

function [traindata_para_cell, testdata_para_cell, Measure_Table, error_rate] = ...
    twoD_NCR_Seg(traindata, trainlabel, testdata, testlabel, SegPoint, Name)
%% ѵ�����Ͳ��Լ���ʼ��
[m1, ~] = size(traindata);             %ѵ������m1��ʱ�����У�ÿ��ʱ��������nά�ȣ�
[m2, n] = size(testdata);              %���Լ���m2��ʱ�����У�ÿ��ʱ��������nά�ȣ�

% disp(['�������ݼ�Ϊ��',Name]);
% disp(['ѵ����������',num2str(m1),'  ���Լ�������',num2str(m2),'  ����ά��Ϊ��',num2str(n)]);

traindata_diff = [zeros(m1, 1) diff(traindata, 1, 2)];     %���ѵ��������һ�ײ�ֺ�����ݣ���ǰ�˲���
testdata_diff = [zeros(m2, 1) diff(testdata, 1, 2)];       %��ֲ��Լ�����һ�ײ�ֺ�����ݣ���ǰ�˲���

%% ʵ������ʱ��ε���ģ�������仯

w = length(SegPoint) + 1; % wΪ�ֶ���
SegPoint = [1 SegPoint n]; % �ֶε������ͷ��β����

traindata_para_cell = cell(m1, 1);    %ѵ������ģ�Ͳ�������ÿ��Ԫ���ں���һ��w*6�ľ��󣬴���һ��ʱ�����ݵ�w�ε���ģ�Ͳ���
testdata_para_cell  = cell(m2, 1);    %ѵ������ģ�Ͳ�������ÿ��Ԫ���ں���һ��w*6�ľ��󣬴���һ��ʱ�����ݵ�w�ε���ģ�Ͳ���

for j = 1: m1 % ��ѵ�����ĵ�j�������Ʋ���
    Cloudparameter = zeros(w, 6);          %��j�����е��Ʋ�����ÿһ����һ��
    for k = 1: w % ���k�ε��Ʋ���
        x1 = traindata(j, SegPoint(k) + (k ~= 1): SegPoint(k + 1));       %x1����ѵ����ʱ��ε��м��������j��ʱ�����еĵ�k��
        x2 = traindata_diff(j, SegPoint(k) + (k ~= 1): SegPoint(k + 1));  %x2������ѵ����ʱ��ε��м����
        [Ex1, En1, He1, Ex2, En2, He2] = backward_twoD_Clouds(x1, x2);   %���ö�ά������ģ���Ӻ���
        Cloudparameter(k, :) = [Ex1, En1, He1, Ex2, En2, He2];
    end
    traindata_para_cell{j} = Cloudparameter;
end

for j = 1: m2 % ����Լ��ĵ�j�������Ʋ���
    Cloudparameter = zeros(w, 6);          %�м����
    for k = 1: w
        x1 = testdata(j, SegPoint(k) + (k ~= 1): SegPoint(k + 1));       %x1����ѵ����ʱ��ε��м����
        x2 =  testdata_diff(j, SegPoint(k) + (k ~= 1): SegPoint(k + 1));  %x2������ѵ����ʱ��ε��м����
        [Ex1, En1, He1, Ex2, En2, He2] = backward_twoD_Clouds(x1, x2);   %���ö�ά������ģ���Ӻ���
        Cloudparameter(k,:) = [Ex1, En1, He1, Ex2, En2, He2];
    end
    testdata_para_cell{j} = Cloudparameter;
end

%% ʵ������ʱ�������ƾ������ 
Measure_Table = zeros(m2, m1); %������󣬵�M�У�N�У������ŵ�M�����Լ�ʱ���������N��ѵ����ʱ�����ݵ��ص������ֵԽ���������ʱ���������ƶ�Խ��

for i=1: m2       %i������Լ�
    for j = 1: m1   %j����ѵ����
        S = 0;     %S�����м����
        for k = 1: w
            C1_parameter = traindata_para_cell{j}(k,:);
            C2_parameter = testdata_para_cell{i}(k,:);
            [TOM] = Measure_for_twoD_Clouds(C1_parameter, C2_parameter);
            S = S + TOM;
        end
        Measure_Table(i, j) = sqrt(S / w) ;
    end
end

%% ʵ��1NN�㷨 
% K = 1; %1NN
[~, I] = sort(Measure_Table,2,'descend'); %ÿ�����ƶ��ɴ�С����
Predictlabel = zeros(m2, 1);              %1NN�㷨��Ԥ���ǩ
for j = 1: m2
    Predictlabel(j, 1) = trainlabel(I(j, 1)) ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Step4  END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%%%%%%%%%%%%%%%%%%%%%%%% Step5��error rate���� %%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    
    
    
    

