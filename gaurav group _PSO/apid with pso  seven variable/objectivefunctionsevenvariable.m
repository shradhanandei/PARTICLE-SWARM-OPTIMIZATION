



function [iae]=objectivefunctionsevenvariable(p)



h=0.1; 
t = 0:h:16;  
tf=160 ;


w1=1;
w2=0;

delay=0.2;
de=2;



for j=1:10


        kp=p(j,1);
        ki=p(j,2);
        kd=p(j,3);

        k1=p(j,4);
        k2=p(j,5);
        k3=p(j,6);
        k4=p(j,7);


y = zeros(1,length(t)+100); 
x = zeros(1,length(t)+100);

u = zeros(1,length(t)+100);


e = zeros(1,length(t)+100);

kpp = zeros(1,length(t)+100); 
kii = zeros(1,length(t)+100);
kdd = zeros(1,length(t)+100);


v = zeros(1,length(t)+100);
g = zeros(1,length(t)+100);


r= 1;

setpoint=r;
y(1) = 0;
x(1)=0 ;
e(1)=r - y(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



u(1)=kp*e(1)+ki*sum(e);

kpp(1)=0;
kii(1)=0;
kdd(1)=0;
v(1)=0;


F_xy = @(x) -2*x; 

for i = 1:de


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

kii(i+1)=ki*(k4+k2*v(i+1));

kdd(i+1)=kd*(1+k3*abs(v(i+1)));


u(i+1)=kpp(i+1)*e(i+1)+ kii(i+1)*sum(e)+ kdd(i+1)*ed;

end




F_xy = @(u,y,x) u-y-2*x;


for i = de+1:(tf/2)



k_1 = F_xy(u(i-de),y(i),x(i));

k_2 = F_xy(u(i-de)+0.5*h,y(i)+0.5*h,x(i)+0.5*h*k_1);

k_3 = F_xy((u(i-de)+0.5*h),(y(i)+0.5*h),(x(i)+0.5*h*k_2));

k_4 = F_xy((u(i-de)+h ) ,(y(i)+h) ,(x(i)+k_3*h));

x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;

y(i+1)=y(i)+h*x(i+1); 


e(i+1)=r-y(i+1);

ed=e(i+1)-e(i);  

er=sum(e);

v(i+1)=(e(i+1)/r)*(ed/r);

kpp(i+1)=kp*(1+k1*abs(v(i+1)));

kii(i+1)=ki*(k4+k2*v(i+1));

kdd(i+1)=kd*(1+k3*abs(v(i+1)));


u(i+1)=kpp(i+1)*e(i+1)+ kii(i+1)*sum(e)+ kdd(i+1)*ed;




end 


u((tf/2)+1)= -20;

F_xy = @(u,y,x) u-y-2*x;


for i = (tf/2)+1:length(t)


k_1 = F_xy(u(i-de),y(i),x(i));

k_2 = F_xy(u(i-de)+0.5*h,y(i)+0.5*h,x(i)+0.5*h*k_1);

k_3 = F_xy((u(i-de)+0.5*h),(y(i)+0.5*h),(x(i)+0.5*h*k_2));

k_4 = F_xy((u(i-de)+h ) ,(y(i)+h) ,(x(i)+k_3*h));

x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;

y(i+1)=y(i)+h*x(i+1); 


e(i+1)=r-y(i+1);

ed=e(i+1)-e(i);

er=sum(e);

v(i+1)=(e(i+1)/r)*(ed/r);

kpp(i+1)=kp*(1+k1*abs(v(i+1)));

kii(i+1)=ki*(k4+k2*v(i+1));

kdd(i+1)=kd*(1+k3*abs(v(i+1)));


u(i+1)=kpp(i+1)*e(i+1)+ kii(i+1)*sum(e)+ kdd(i+1)*ed;

end 




z=y(1:(tf+1));


figure(j+1)

plot(t,z) 


xlabel('Time t ')

ylabel('Response y')

title(' exp(-0.2s) /(s+1)^2  loop ')

grid on 


%%%%%%%%   IAE   AND ITAE    %%%%%

iaeiae(j,1)=0.1*sum(abs(e));


for i= 1:length(t)

g(i)=0.01*i*e(i);


end

itaeitae(j,1)=sum(abs(g));

end


iae = w1*iaeiae+w2*itaeitae;



