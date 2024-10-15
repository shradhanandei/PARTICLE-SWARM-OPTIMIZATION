%%tf=exp(-0.2s)/s(s+1)
%%ultimate gain and ultimate time period calulation 

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
r=1
y(1) = 0;
x(1)=0 ;
e(1)=r- y(1);

%%%%%%%%%%%%

kp=5.09;

u(1)=kp*e(1);
 
delay=0.2
de=2



u(1)=kp*(e(1))

F_xy = @(x) -x; 


for i = 1:2
  
 
   k_1 = F_xy(x(i));
   
   k_2 = F_xy(x(i)+0.5*h*k_1);
   
   k_3 = F_xy((x(i)+0.5*h*k_2));
   
   k_4 = F_xy((x(i)+k_3*h));

  
    x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    
   y(i+1)=y(i)+h*x(i+1);

e(i+1)=1-y(i+1);

u(i+1)=kp*e(i+1);

end



F_xy = @(u,x) u-x;


for i = 3:length(t)
    
    k_1 = F_xy(u(i-2),x(i));
    
    
    k_2 = F_xy(u(i-2)+0.5*h,x(i)+0.5*h*k_1);
    
    
    k_3 = F_xy((u(i-2)+0.5*h),(x(i)+0.5*h*k_2));
    
    
    k_4 = F_xy((u(i-2)+h ) ,(x(i)+k_3*h)); 
    
    
    x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    
   y(i+1)=y(i)+h*x(i+1);
   
   e(i+1)=1-y(i+1);
  u(i+1)=kp*e(i+1);   
end 

z=y(1:301)
plot(t,z) 
xlabel('value of time t')
ylabel('value of y')
title(' calculation of ultimate gain,k =5.09 or(5.11 approx) ultimate period=2.9sec, second order integrating system tf=exp(-0.2s)/s(s+1) ')
grid on