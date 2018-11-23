function [OD,u]=Measure_for_oneD_Clouds(Ex1,En1,He1,Ex2,En2,He2)
    
if En1 == 0
    En1 = 10^-30 ;
end

if En2 == 0 
    En2 = 10^-30 ;
end

SupC1=Ex1+3*En1;
InfC1=Ex1-3*En1;
SupC2=Ex2+3*En2;
InfC2=Ex2-3*En2;

OD = 2*(min(SupC1,SupC2)-max(InfC1,InfC2))/(6*En1+6*En2) ;

if  En1==En2
    x1=(Ex1*En2+Ex2*En1)/(En1+En2);
    n=1; %����ֻ��һ�����㣬��En1=En2ʱ���Ž���һ����
else
    a=(Ex2*En1-Ex1*En2)/(En1-En2);
    b=(Ex1*En2+Ex2*En1)/(En1+En2);
    x1=min(a,b);
    x2=max(a,b);   
    n=2; %������2������
end
%ֻ��һ������ʱ��uֵ
e = 10^-20; %������ȡ���Ⱥŵ��ж�
if  n==1
    if  (max(InfC1,InfC2)-e) < x1 && x1< (min(SupC1,SupC2)+e)
        u = exp( - (x1 - Ex1) ^2 / (2* En1^2) );
    else
        u=0;
    end
end

%����������ʱ��uֵ
if  n==2
    u1 = 0; u2 = 0; %��ʼ�����������uֵ
    if  (max(InfC1,InfC2)-e) < x1 && x1< (min(SupC1,SupC2)+e)
        u1 = exp( - (x1 - Ex1) ^2 / (2* En1^2) );    
    end
    if  (max(InfC1,InfC2)-e) < x2 && x2< (min(SupC1,SupC2)+e)
        u2 = exp( - (x2 - Ex2) ^2 / (2* En2^2) );    
    end
    u = max(u1,u2); 
end
end