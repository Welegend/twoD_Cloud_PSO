%{ 
调用函数：Measure_for_oneD_Clouds.m;
%}

function [TOM]=Measure_for_twoD_Clouds(C1_parameter,C2_parameter)
   %C1_parameter=[Ex1,En1,He1,Ex2,En2,He2] 代表一个二维云模型
   %C2_parameter=[Ex3,En3,He3,Ex4,En5,He6] 代表另一个二维云模型
%一朵二维云的6个参数   
Ex1=C1_parameter(1); En1=C1_parameter(2); He1=C1_parameter(3);
Ex2=C1_parameter(4); En2=C1_parameter(5); He2=C1_parameter(6);
%另一朵二维云的6个参数   
Ex3=C2_parameter(1); En3=C2_parameter(2); He3=C2_parameter(3);
Ex4=C2_parameter(4); En4=C2_parameter(5); He4=C2_parameter(6);

[OD1,u1]=Measure_for_oneD_Clouds(Ex1,En1,He1,Ex3,En3,He3);
[OD2,u2]=Measure_for_oneD_Clouds(Ex2,En2,He2,Ex4,En4,He4);
OM1=OD1*u1;
OM2=OD2*u2;
TOM=OM1*OM2;
%TOM=OM1+OM2;
end