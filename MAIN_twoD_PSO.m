% �������ܣ���PSO�㷨����ֶ���2-9ʱ����ѷֶε㣬һ��ֻ�ܲ������е�һ�����ݼ�
% ���룺���ݼ�����Ԫ������
% �������Ӧ�ֶ���w�£���ѷֶε�SegPoint_cell����С�����error_rate��������w������ʵ�������

function [w, SegPoint_cell, error_rate] = MAIN_twoD_PSO(newName_Dataset)
%% �������ݼ�
i = 1; % ����newName_Dataset�еĵ�һ�����ݼ�
Name = newName_Dataset{i, 2}; % ���ݼ�����
Path = 'F:\�о���\�����ھ������\ʵ��ģ��\��ά��ģ��\UCR_TS_Archive_2015';
eval(['load ', Path, '\', newName_Dataset{i, 2}, '\',Name, '_TRAIN']);
eval(['load ', Path, '\', newName_Dataset{i, 2}, '\',Name, '_TEST']);

%% �����ݼ��õ�w_PSO���������
Rowtraindata = eval([newName_Dataset{i, 2}, '_TRAIN']);
Rowtestdata = eval([newName_Dataset{i, 2}, '_TEST']);

Rowtraindata = sortrows(Rowtraindata, 1) ;             %�������Ų���һ�𣨰���һ�ж��н��������У�
Rowtestdata = sortrows(Rowtestdata, 1) ;

traindata = Rowtraindata(:, 2: end);       %ѵ��������ÿ�д���һ��ʱ������
testdata = Rowtestdata(:, 2: end);         %���Լ�����ÿ�д���һ��ʱ������
trainlabel = Rowtraindata(:, 1);          %ѵ������ǩ
testlabel = Rowtestdata(:, 1);            %���Լ���ǩ

%% ��PSO�㷨��ֶ�����2��9����ѷֶε�������
w = ([4, 8])'; % wΪҪѰ�ŵķֶ���
SegPoint_cell = cell(length(w), 1);
error_rate = zeros(length(w), 1);
for i = 1: length(w)
    [SegPoint_cell{i}, error_rate(i)] = w_PSO(traindata, trainlabel, testdata, testlabel, w(i), Name);
end

end




