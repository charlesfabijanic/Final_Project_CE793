% Code Thwaites method
%================================%
%========= Setup ================%
clc;clear all;
fun=@(x) (2*sin(x)).^5;
mu=10^-5;
Re_inf=10^5;
a=1;
U_inf=1;
x_max=105*pi/180;
x=[0.01:.001:x_max];
theta=zeros(1,length(x));
lambda=zeros(1,length(x));
S_lambda=zeros(1,length(x));
Cf=zeros(1,length(x));
I=zeros(1,length(x));
dI=0;
xmin=0;
for i=1:length(x)
    I(i)=integral(fun,0,x(i));
    theta(i)=sqrt(0.45*mu/(f(x(i)))^6*I(i));
    lambda(i)=theta(i)^2*potential_diff(x(i))/mu;
    S_lambda(i)=(lambda(i)+0.09)^0.62;
    Cf(i)=2*mu*S_lambda(i)/(f(x(i))*theta(i));
end
ind=find(lambda<=-0.089);
x(ind(1))*180/pi
figure (1)
plot(x*180/pi,theta,'g','LineWidth',2)
xlabel('\beta','FontSize',11,'FontWeight','bold');
ylabel('\theta','FontSize',11,'FontWeight','bold');
title(' \theta Vs \beta');

figure(2)
plot(x*180/pi,lambda,'m','LineWidth',2)
xlabel('\beta','FontSize',11,'FontWeight','bold');
ylabel('\lambda','FontSize',11,'FontWeight','bold');
title(' \lambda Vs \beta ')

figure(3)
plot(x*180/pi,Cf,'k','LineWidth',2)
xlabel('\beta','FontSize',11,'FontWeight','bold');
ylabel('Cf','FontSize',11,'FontWeight','bold');
title(' Cf Vs \beta ')
