
%%%% calculation of ultimate gain and time period 
%%%%tf=exp(-0.2s)/(s2+2s+1)

clc;                                              
clear all;

h=0.1;                                           
t = 0:h:30; 
y = zeros(1,length(t)); 
u = zeros(1,length(t));
e = zeros(1,length(t)); 
x = zeros(1,length(t));
i = zeros(1,length(t));
input = 1;
y(1) = 0;
x(1)=0 ;
e(1)=input - y(1);



kp=10.5;

u(1)=kp*e(1);

F_xy = @(x) -2*x; 

for i = 1:2
     
%y and u taken as time input and x as output

k_1 = F_xy(x(i));
k_2 = F_xy(x(i)+0.5*h*k_1);
k_3 = F_xy((x(i)+0.5*h*k_2));
k_4 = F_xy((x(i)+k_3*h)); 


x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;

y(i+1)=y(i)+h*x(i+1);

e(i+1)=1-y(i+1);

u(i+1)=kp*e(i+1);

end

F_xy = @(u,y,x) u-y-2*x;


for i = 3:length(t)

k_1 = F_xy(u(i-2),y(i),x(i));

k_2 = F_xy(u(i-2)+0.5*h,y(i)+0.5*h,x(i)+0.5*h*k_1);

k_3 = F_xy((u(i-2)+0.5*h),(y(i)+0.5*h),(x(i)+0.5*h*k_2));

k_4 = F_xy((u(i-2)+h ) ,(y(i)+h) ,(x(i)+k_3*h));


x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;

y(i+1)=y(i)+h*x(i+1);

e(i+1)=1-y(i+1);

u(i+1)=kp*e(i+1);

end  


z=y(1:301)
plot(t,z) 
xlabel('value of time t')
ylabel('value of y')
title(' calculation of ultimate gain,k =10.5 ultimate period=2.033sec, second order system tf=exp(-0.2s)/(s2+2s+1) ')
grid on




