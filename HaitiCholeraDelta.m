
%%%%%Sensitivity analysis of Simulation 1 of Haiti Cholera paper 5/19/2020
%%%sensitivity analysis to beta
%%%%%
%%%%%%%%
P=5;%%%P is the percentage change of beta
DT=0.01/2; %% Time interval lenght=DT day
T=100;  %%Max number of days
Tmax=T/DT;    %%%Time iterations
DX=1;  %%%Space interval length=1 km
Xmax=1000; %%%Space length= Xmax*Dx


%%%%%%%
%%Initial conditions and variables

%%Susceptible
S1=zeros(Tmax,1);S2=zeros(Tmax,1);S3=zeros(Tmax,1);
S1(1)=10^5; S2(1)=10^5;S3(1)=10^5;
%%Infected
I1=zeros(Tmax,1);I2=zeros(Tmax,1);I3=zeros(Tmax,1);
I1(1)=0; I2(1)=0; I3(3)=0;
%%Bacteria in community
B1=zeros(Tmax,1);B2=zeros(Tmax,1);B3=zeros(Tmax,1);
B1(1)=0; B2(1)=0; B3(3)=0;
%%Bacteria in water
b=zeros(Xmax, Tmax);

%%%%%%%%%
%%Parameters
%beta=0.025;
beta=1;%beta=0.03;
delta=0.23;
gamma=0.19; 
r=0.18*10^-5;

mu=5*10^-2; tildemu=mu*100;

 rho=0.5; tilderho=0.5*10^-2;

 epsilon=0.1;
d=1;
q=50;

b(495, 1)=10^-1; 

for m=1:(Tmax-1)
    S1(m+1)=S1(m)-DT*beta*S1(m)*B1(m)/(1+B1(m));
    S2(m+1)=S2(m)-DT*beta*S2(m)*B2(m)/(1+B2(m));
    S3(m+1)=S3(m)-DT*beta*S3(m)*B3(m)/(1+B3(m));
    
    I1(m+1)=I1(m)+DT*(beta*S1(m)*B1(m)/(1+B1(m))-gamma*I1(m));
    I2(m+1)=I2(m)+DT*(beta*S2(m)*B2(m)/(1+B2(m))-gamma*I2(m));
    I3(m+1)=I3(m)+DT*(beta*S3(m)*B3(m)/(1+B3(m))-gamma*I3(m));
    
    B1(m+1)=B1(m)+DT*r*I1(m)+DT*tildemu*0.5*(b(500,m)+b(501, m))-delta*DT*B1(m)-rho*DT*B1(m);  %%%%First site at x=500
    B2(m+1)=B2(m)+DT*r*I2(m)+DT*tildemu*0.5*(b(600,m)+b(601, m))-delta*DT*B2(m)-rho*DT*B2(m);  %%%%Second site at x=600
    B3(m+1)=B3(m)+DT*r*I3(m)+DT*tildemu*0.5*(b(620,m)+b(621, m))-delta*DT*B3(m)-rho*DT*B3(m);  %%%%Third site at x=620
    
    for n=2: (Xmax-1)
        if n==500||n==501
            b(n, m+1)=b(n,m)+epsilon/(DX)^2*DT*(b(n+1, m)-2*b(n, m)+b(n-1, m))-q/DX*DT*(b(n, m)-b(n-1, m))-delta*DT*b(n,m)+DT*tilderho/d*B1(m)-DT*mu*b(n, m);
        elseif n==600||n==601 
            b(n, m+1)=b(n,m)+epsilon/(DX)^2*DT*(b(n+1, m)-2*b(n, m)+b(n-1, m))-q/DX*DT*(b(n, m)-b(n-1, m))-delta*DT*b(n,m)+DT*tilderho/d*B2(m)-DT*mu*b(n, m);
        elseif n==620||n==621 
            b(n, m+1)=b(n,m)+epsilon/(DX)^2*DT*(b(n+1, m)-2*b(n, m)+b(n-1, m))-q/DX*DT*(b(n, m)-b(n-1, m))-delta*DT*b(n,m)+DT*tilderho/d*B3(m)-DT*mu*b(n, m);
        else
            b(n, m+1)=b(n,m)+epsilon/(DX)^2*DT*(b(n+1, m)-2*b(n, m)+b(n-1, m))-q/DX*DT*(b(n, m)-b(n-1, m))-delta*DT*b(n,m);
        end
    end
end

% figure
% plot(0:DT:(T-DT), I1)
% hold on 
% plot(0:DT:(T-DT), I2)
% plot(0:DT:(T-DT), I3)
% hold off 

% h=figure
% plot(0:DT:(T-DT), S1, 'LineWidth', 1)
% hold on
% xlabel('days')
% ylabel('cases')
% plot(0:DT:(T-DT), S2, 'LineWidth', 1) 
% plot(0:DT:(T-DT), S3, 'LineWidth', 1) 
% legend('Community 1', 'Community 2', 'Community 3') 
% saveas(h, 'Ex1Susceptible', 'epsc')
% hold off

h=figure
%subplot(1,3, 1)
case1=beta*S1.*B1./(1+B1);  %%Onset new cases in 1st community
plot(0:DT:(T-DT), case1, 'LineWidth', 1) 
xlabel('days')
ylabel('cases')
%title('Community 1')

hold on
%subplot(1,3, 2)
case2=beta*S2.*B2./(1+B2);  %%Onset new cases in 1st community
plot(0:DT:(T-DT), case2, 'LineWidth', 1) 
%xlabel('days')
%ylabel('cases')

case3=beta*S3.*B3./(1+B3);  %%Onset new cases in 1st community
plot(0:DT:(T-DT), case3, 'LineWidth', 1) 

%title('Community 2')

% legend('Community 1', 'Community 2', 'Community 3') 
% saveas(h, 'Ex1Reported', 'epsc')


T1=S1(1)+S2(1)+S3(1)-S1(Tmax)-S2(Tmax)-S3(Tmax);



%%%%%
%%Space Time
%%%%%%%%
DT=0.01/2; %% Time interval lenght=DT day
T=100;  %%Max number of days
Tmax=T/DT;    %%%Time iterations
DX=1;  %%%Space interval length=1 km
Xmax=1000; %%%Space length= Xmax*Dx


%%%%%%%
%%Initial conditions and variables

%%Susceptible
S1=zeros(Tmax,1);S2=zeros(Tmax,1);S3=zeros(Tmax,1);
S1(1)=10^5; S2(1)=10^5;S3(1)=10^5;
%%Infected
I1=zeros(Tmax,1);I2=zeros(Tmax,1);I3=zeros(Tmax,1);
I1(1)=0; I2(1)=0; I3(3)=0;
%%Bacteria in community
B1=zeros(Tmax,1);B2=zeros(Tmax,1);B3=zeros(Tmax,1);
B1(1)=0; B2(1)=0; B3(3)=0;
%%Bacteria in water
b=zeros(Xmax, Tmax);

%%%%%%%%%
%%Parameters
%beta=0.025;
beta=1;%beta=0.03;
delta=0.23*(1+P/100);
gamma=0.19; 
r=0.18*10^-5;

mu=5*10^-2; tildemu=mu*100;

rho=0.5; tilderho=0.5*10^-2;

epsilon=0.1;
d=1;
q=50;

b(495, 1)=0.1;  

for m=1:(Tmax-1)
    S1(m+1)=S1(m)-DT*beta*S1(m)*B1(m)/(1+B1(m));
    S2(m+1)=S2(m)-DT*beta*S2(m)*B2(m)/(1+B2(m));
    S3(m+1)=S3(m)-DT*beta*S3(m)*B3(m)/(1+B3(m));
    
    I1(m+1)=I1(m)+DT*(beta*S1(m)*B1(m)/(1+B1(m))-gamma*I1(m));
    I2(m+1)=I2(m)+DT*(beta*S2(m)*B2(m)/(1+B2(m))-gamma*I2(m));
    I3(m+1)=I3(m)+DT*(beta*S3(m)*B3(m)/(1+B3(m))-gamma*I3(m));
    
    B1(m+1)=B1(m)+DT*r*I1(m)+DT*tildemu*0.5*(b(500,m)+b(501, m))-delta*DT*B1(m)-rho*DT*B1(m);  %%%%First site at x=500
    B2(m+1)=B2(m)+DT*r*I2(m)+DT*tildemu*0.5*(b(600,m)+b(601, m))-delta*DT*B2(m)-rho*DT*B2(m);  %%%%Second site at x=600
    B3(m+1)=B3(m)+DT*r*I3(m)+DT*tildemu*0.5*(b(620,m)+b(621, m))-delta*DT*B3(m)-rho*DT*B3(m);  %%%%Third site at x=620
    
    for n=2: (Xmax-1)
        if n==500||n==501
            b(n, m+1)=b(n,m)+epsilon/(DX)^2*DT*(b(n+1, m)-2*b(n, m)+b(n-1, m))-q/DX*DT*(b(n, m)-b(n-1, m))-delta*DT*b(n,m)+DT*tilderho/d*B1(m)-DT*mu*b(n, m);
        elseif n==600||n==601 
            b(n, m+1)=b(n,m)+epsilon/(DX)^2*DT*(b(n+1, m)-2*b(n, m)+b(n-1, m))-q/DX*DT*(b(n, m)-b(n-1, m))-delta*DT*b(n,m)+DT*tilderho/d*B2(m)-DT*mu*b(n, m);
        elseif n==620||n==621 
            b(n, m+1)=b(n,m)+epsilon/(DX)^2*DT*(b(n+1, m)-2*b(n, m)+b(n-1, m))-q/DX*DT*(b(n, m)-b(n-1, m))-delta*DT*b(n,m)+DT*tilderho/d*B3(m)-DT*mu*b(n, m);
        else
            b(n, m+1)=b(n,m)+epsilon/(DX)^2*DT*(b(n+1, m)-2*b(n, m)+b(n-1, m))-q/DX*DT*(b(n, m)-b(n-1, m))-delta*DT*b(n,m);
        end
    end
end

% figure
% plot(0:DT:(T-DT), I1)
% hold on 
% plot(0:DT:(T-DT), I2)
% plot(0:DT:(T-DT), I3)
% % hold off 
% 
% h=figure
% plot(0:DT:(T-DT), S1, 'LineWidth', 1)
% hold on
% xlabel('days')
% ylabel('cases')
% plot(0:DT:(T-DT), S2, 'LineWidth', 1) 
% plot(0:DT:(T-DT), S3, 'LineWidth', 1) 
% legend('Community 1', 'Community 2', 'Community 3') 
% saveas(h, 'Ex1Susceptible', 'epsc')
% hold off
% 
% h=figure
%subplot(1,3, 1)
case1=beta*S1.*B1./(1+B1);  %%Onset new cases in 1st community
plot(0:DT:(T-DT), case1,'--', 'LineWidth', 1) 
xlabel('days')
ylabel('cases')
%title('Community 1')


%subplot(1,3, 2)
case2=beta*S2.*B2./(1+B2);  %%Onset new cases in 1st community
plot(0:DT:(T-DT), case2, '--','LineWidth', 1) 
%xlabel('days')
%ylabel('cases')

case3=beta*S3.*B3./(1+B3);  %%Onset new cases in 1st community
plot(0:DT:(T-DT), case3, '--','LineWidth', 1) 

%title('Community 2')

% legend('Community 1, \beta_j=1', 'Community 2, \beta_j=1', 'Community 3, \beta_j=1', 'Community 1, \beta_j=1.2', 'Community 2, \beta_j=1.2', 'Community 3, \beta_j=1.2') 
% saveas(h, 'sens_beta', 'epsc')

T2=S1(1)+S2(1)+S3(1)-S1(Tmax)-S2(Tmax)-S3(Tmax);
output=(T2-T1)/T1/(P/100)


