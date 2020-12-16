clc;clear all;close all;
nu=10^-5;
angle_step=0.5;
dx=angle_step*pi/180; % Step size in x direction(radial)
dy=.0001; %Step size in y direction
x=0:angle_step:120; % Discretisizing x
y=0:dy:.015; % Discretisizing y
M=length(x);
N=length(y);
U=potential(x*pi./180); % Velocity Distribution in boundary layer
u=zeros(N,M);
v=zeros(N,M);
% Boundary Condition %
u(1,:)=0;v(1,:)=0;
u(2:end-1,1)=.001;
v(2:end-1,1)=-.0001;
u(end,:)=U;  
u(end,1)=.001;
P=zeros(N-1,M);
Q=zeros(N-1,M);
C=zeros(N-1,M);
alpha=zeros(N-1,M);
beta=zeros(N-1,M);
theta=zeros(1,M);
%d_theta=zeros(N,M-1);
Cf=zeros(1,M);
for m=1:M-1
    alpha(2,m)=nu*dx/(u(2,m)*dy^2);
    beta(2,m)=v(2,m)*dx/(2*u(2,m)*dy);
    C(2,m)=u(2,m)-beta(2,m)*(u(3,m)-u(1,m))+(u(N,m+1)^2-u(N,m)^2)/...
        (2*u(2,m));
    P(2,m)=alpha(2,m)/(1+2*alpha(2,m));
    Q(2,m)=C(2,m)/(1+2*alpha(2,m));
    for n=3:(N-1)
        alpha(n,m)=nu*dx/(u(n,m)*dy^2);
        beta(n,m)=v(n,m)*dx/(2*u(n,m)*dy);
        C(n,m)=u(n,m)-beta(n,m)*(u(n+1,m)-u(n-1,m))+(u(N,m+1)^2-u(N,m)^2)/...
            (2*u(n,m));
        P(n,m)=alpha(n,m)/(1+2*alpha(n,m)-alpha(n,m)*P(n-1,m));
        Q(n,m)=(C(n,m)+alpha(n,m)*Q(n-1,m))/(1+2*alpha(n,m)-...
            alpha(n,m)*P(n-1,m));
    end
    u(N-1,m+1)=(alpha(N-1,m)*u(N,m+1)+alpha(N-1,m)*Q(N-2,m)+C(N-1,...
        m))/(1+2*alpha(N-1,m)-alpha(N-1,m)*P(N-2,m));
    for n=N-2:-1:2
        u(n,m+1)=P(n,m)*u(n+1,m+1)+Q(n,m);
    end
    for n=2:N
        v(n,m+1)=v(n-1,m+1)-(dy*(u(n,m+1)-u(n,m)+u(n-1,m+1)-u(n-1,m)))/(2*dx);
    end
    for n=1:N
    d_theta(n,m)=u(n,m+1)/u(N,m+1)*(1-u(n,m+1)/u(N,m+1));
    end
    
    theta=trapz(dy,d_theta);
        Cf(m)=2*nu*(u(2,m)-u(1,m))/dy;
        if Cf(m)<0
            break
        end
        
end

    fprintf('The separation is Predicted at %4.2f',x(m+1))
    figure(1)
    plot(x(:,1:m-1),theta(1,1:m-1),'LineWidth',2)
    xlabel('\beta','FontSize',12,'FontWeight','bold')
    ylabel('\theta','FontSize',12,'FontWeight','bold')
    title(' \theta Vs \beta')
   
    figure(2)
    plot(x(1,1:m-1),Cf(1,1:m-1),'LineWidth',2)
    xlabel('\beta','FontSize',12,'FontWeight','bold')
    ylabel('Cf','FontSize',12,'FontWeight','bold')
    title(' Cf Vs \beta')
    
    figure(3)
    plot(u(:,41),y,'-.r','LineWidth',2)
    hold on
    plot(u(:,121),y,'--b','LineWidth',2)
    hold on
    plot(u(:,181),y,':g','LineWidth',2)
    hold on
    plot(u(:,m),y,'m','LineWidth',2)
    hold on
    legend('\beta=20^{0}','\beta=60^{0}','\beta=90^{0}','\beta=107^{0}'...
        ,'Location','northwest')
    
    xlabel('u/U_{0}','FontSize',12,'FontWeight','bold')
    ylabel('y/R','FontSize',12,'FontWeight','bold')
    title('Velocity Profile')
