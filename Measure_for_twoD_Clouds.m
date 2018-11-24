%{ 
���ú�����Measure_for_oneD_Clouds.m;
���룺6��ѵ����������6�����Լ����������Ǿ������������ݼ����и�������������������
���������TOM
%}

function TOM = Measure_for_twoD_Clouds( ...
    train_Ex1, train_En1, train_He1, train_Ex2, train_En2, train_He2, ...
    test_Ex1, test_En1, test_He1, test_Ex2, test_En2, test_He2)

% ע�����������ѵ�����Ͳ��Լ����佻���
[OD1, u1] = Measure_for_oneD_Clouds(train_Ex1, train_En1, train_He1, test_Ex1, test_En1, test_He1);
[OD2, u2] = Measure_for_oneD_Clouds(train_Ex2, train_En2, train_He2, test_Ex2, test_En2, test_He2);

% �������
OM1 = OD1 .* u1;
OM2 = OD2 .* u2;
TOM = OM1 .* OM2;

end