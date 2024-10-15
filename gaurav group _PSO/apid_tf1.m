
%%%%APID  adaptive PID CONTROLLER  Response of second order system TF=exp(-0.2s)/s2+2s+1

clc;                                               
clear all;

h=0.1;  
t = 0:h:16;   
tf=16/h;

y = zeros(1,length(t)); 
u = zeros(1,length(t));
e = zeros(1,length(t)); 
x = zeros(1,length(t));

kpp = zeros(1,length(t)); 
kii = zeros(1,length(t));
kdd = zeros(1,length(t));

v = zeros(1,length(t));
x = zeros(1,length(t));

r= 1;

setpoint=r;
y(1) = 0;
x(1)=0 ;
e(1)=r - y(1);


ku=10.5;
tu=2.033;

kp=0.6*ku
ti=tu/2
td=tu/8

ki=kp*(0.1/ti)
kd=kp*(td/0.1)

k1=1
k2=1
k3=12

u(1)=kp*e(1)+ki*sum(e);

 kpp(1)=0;
 kii(1)=0;
 kdd(1)=0;
 v(1)=0;


F_xy = @(x) -2*x; 

for i = 1:2


k_1 = F_xy(x(i));
k_2 = F_xy(x(i)+0.5*h*k_1);
k_3 = F_xy((x(i)+0.5*h*k_2));
k_4 = F_xy((x(i)+k_3*h));

x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;

y(i+1)=y(i)+h*x(i+1);

e(i+1)=r-y(i+1);

er=sum(e);

ed=e(i+1)-e(i);

v(i+1)=(e(i+1)/r)*(ed/r);

kpp(i+1)=kp*(1+k1*abs(v(i+1)));
kii(i+1)=ki*(0.3+k2*v(i+1));
kdd(i+1)=kd*(1+k3*abs(v(i+1)));


u(i+1)=kpp(i+1)*e(i+1)+ kii(i+1)*sum(e)+ kdd(i+1)*ed;

end


F_xy = @(u,y,x) u-y-2*x;

for i = 3:(tf/2)
k_1 = F_xy(u(i-2),y(i),x(i));

k_2 = F_xy(u(i-2)+0.5*h,y(i)+0.5*h,x(i)+0.5*h*k_1);

k_3 = F_xy((u(i-2)+0.5*h),(y(i)+0.5*h),(x(i)+0.5*h*k_2));

k_4 = F_xy((u(i-2)+h ) ,(y(i)+h) ,(x(i)+k_3*h));

x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;

y(i+1)=y(i)+h*x(i+1); 

e(i+1)=r-y(i+1);
ed=e(i+1)-e(i);  
er=sum(e);
v(i+1)=(e(i+1)/r)*(ed/r);
kpp(i+1)=kp*(1+k1*abs(v(i+1)));

kii(i+1)=ki*(0.3+k2*v(i+1));

kdd(i+1)=kd*(1+k3*abs(v(i+1)));


u(i+1)=kpp(i+1)*e(i+1)+ kii(i+1)*sum(e)+ kdd(i+1)*ed;

end 


u((tf/2)+1)= -14;

F_xy = @(u,y,x) u-y-2*x;


for i = (tf/2)+1:length(t)
    k_1 = F_xy(u(i-2),y(i),x(i));
    
    k_2 = F_xy(u(i-2)+0.5*h,y(i)+0.5*h,x(i)+0.5*h*k_1);
    
    k_3 = F_xy((u(i-2)+0.5*h),(y(i)+0.5*h),(x(i)+0.5*h*k_2));
    
    k_4 = F_xy((u(i-2)+h ) ,(y(i)+h) ,(x(i)+k_3*h));
    
    x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    
   y(i+1)=y(i)+h*x(i+1); 

e(i+1)=r-y(i+1);
ed=e(i+1)-e(i);  
er=sum(e);
v(i+1)=(e(i+1)/r)*(ed/r);
kpp(i+1)=kp*(1+k1*abs(v(i+1)));

kii(i+1)=ki*(0.3+k2*v(i+1));

kdd(i+1)=kd*(1+k3*abs(v(i+1)));


u(i+1)=kpp(i+1)*e(i+1)+ kii(i+1)*sum(e)+ kdd(i+1)*ed;
    
end 


z=y(1:(tf+1));

plot(t,z) 
hold on 
xlabel('Time t ')
ylabel('Response y')
title('APID  Response of second order system TF=exp(-0.2s)/s2+2s+1) ')
grid on 








   

