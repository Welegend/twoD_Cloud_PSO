function [Ex1,En1,He1,Ex2,En2,He2]=backward_twoD_Clouds(x1,x2)
 %输入两组等长时间序列段
 n=length(x1);
 Ex1=mean(x1);
 Ex2=mean(x2);
 En1=mean(abs(x1-Ex1)) * sqrt(pi/2);
 En2=mean(abs(x2-Ex2)) * sqrt(pi/2);
 He1=sqrt(var(x1)-En1^2);
 He2=sqrt(var(x2)-En2^2);
end