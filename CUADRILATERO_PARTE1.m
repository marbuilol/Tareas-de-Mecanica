close all; clc; close all;
format long
%% Vector phi_2
phi_2=[100
90
80
70
60
50
40
30
20
10
0]';
phi_2=(pi()/180).*phi_2; % phi_2 en radiaes para poder trabajar con ellos
%% Declaración de datos y vectores 

w_21=1; %velocidad angular de la barra 2 que es cte
a=0.1; %unidad de longitud 
phi0=pi()/2; %ángulo inicial de la barra 2
l_1=2*a;
l_2=a;
l_3=a;
l_4=sqrt(2).*a;

N=length(phi_2); %longitud del vector phi_2

t=zeros(1,N); %creación del vector fila del tiempo de dim N

phi_3=zeros(1,N);%creación del vector fila phi_3 de dim N
phi_4=zeros(1,N);%creación del vector fila phi_4 de dim N

w_31=zeros(1,N);%creación del vector fila omega_31 de dim N
w_41=zeros(1,N);%creación del vector fila omega_41 de dim N

dw_31dt=zeros(1,N);%creación del vector fila alpha_3 de dim N
dw_41dt=zeros(1,N);%creación del vector fila alpha_4 de dim N

vB=zeros(1,N); %creación del vector fila vB_3 de dim N
vC=zeros(1,N); %creación del vector fila vC_3 de dim N

aB=zeros(1,N);%creación del vector fila aB_3 de dim N
aC=zeros(1,N);%creación del vector fila aC_3 de dim N
%% Cálculo del tiempo
%se ha tenido en cuenta que el sistema parte en phi_2 cuando es 90º
for i=1:1:N %bucle que calcula los tiempos que tarda en llegar la barra B
            %a los ángulos correspondientes phi_2 estando incialmente en
            %el ángulo phi0
    if i==1
    t(i)=abs((phi0-phi_2(i))/w_21);
    else
    t(i)=(i).*t(1);
    end    
end
%% Cálculo de los ángulos phi_3 y phi_4

%Utilizando las fórmulas adjuntadas en el pdf
for i=1:1:N
A=-l_1+l_2.*cos(phi_2(i));

B=l_2.*sin(phi_2(i));

C=((l_4.^2)-(l_3.^2)-((A.^2)+(B.^2)))./(2.*l_3);

D=((l_4.^2)-(l_3.^2)+((A.^2)+(B.^2)))./(2.*l_4);

phi_3(i)=2.*atan((B-sqrt((A.^2)+(B.^2)-(C.^2)))./(C+A));
phi_4(i)=2.*atan((B-sqrt((A.^2)+(B.^2)-(D.^2)))./(D+A));
end
%paso de radianes a grados sexagesimales
phi_3=(180/pi)*phi_3;
phi_4=(180/pi)*phi_4;

%escogiendo los ángulos menores positivos con respecto al eje X
for i=1:N 
    if phi_3(i)>90
        phi_3(i)=abs(phi_3(i)-180);
    elseif phi_3(i)<0
        phi_3(i)=phi_3(i)+180;
    end
     if phi_3(i)>90
        phi_3(i)=abs(phi_3(i)-180);
    end
end
for i=1:N
    if phi_4(i)>90
        phi_4(i)=abs(phi_4(i)-180);
    elseif phi_4(i)<0
        phi_4(i)=phi_4(i)+180;
    end
    if phi_4(i)>90
        phi_4(i)=abs(phi_4(i)-180);
    end
    
end

%Paso de sexagesimales a radianes
phi_3=(pi()/180)*phi_3;
phi_4=(pi()/180)*phi_4;
%% Cálculo de ángulos de importancia

%ángulo máximo en el que las tres barras pueden girar en sentido antihorario
%instante donde las barras 3 y 4 están alineadas

alpha_1=acos((((1+sqrt(2))^2)-5)/-4); %aquí se comprueba que los 100º
                                      %pertenecen al régimen en el que las
                                      %barras están unidas

%ángulo en el cual la barra 4 cambia de sentido cuando giran horariamente
%hasta 0º
alpha_2=acos(6/8);%en este instante las barras 2 y 3 están alineadas
%% Cálculo de las velocidades angulares y lineales y la aceleracióm de B
%Apoyándome en el planteamiento de las ecuaiones que adjunté
for i=1:1:N
    
    if  phi_2(i)<alpha_1 && phi_2(i)>=(pi()/2)
    a=-w_21*l_2*sin(phi_2(i));
    b=l_3*sin(phi_3(i));
    c=-l_4*sin(phi_4(i));
    d=-w_21*l_2*cos(phi_2(i));
    e=l_3*cos(phi_3(i));
    f=-l_4*cos(phi_4(i));
    w_41(i)=((a*e)-(b*d))/((c*e)-(b*f));
    w_31(i)=((a*f)-(c*d))/((b*f)-(c*e));
    
    elseif phi_2(i)<(pi()/2) && phi_2(i)>=alpha_2
    a=w_21*l_2*sin(phi_2(i));
    b=l_3*sin(phi_3(i));
    c=l_4*sin(phi_4(i));
    d=-w_21*l_2*cos(phi_2(i));
    e=-l_3*cos(phi_3(i));
    f=l_4*cos(phi_4(i));
    w_41(i)=((a*e)-(b*d))/((c*e)-(b*f));
    w_31(i)=((a*f)-(c*d))/((b*f)-(c*e));
    
    elseif phi_2(i)<alpha_2 && phi_2(i)>=0
    a=w_21*l_2*sin(phi_2(i));
    b=l_3*sin(phi_3(i));
    c=-l_4*sin(phi_4(i));
    d=-w_21*l_2*cos(phi_2(i));
    e=-l_3*cos(phi_3(i));
    f=-l_4*cos(phi_4(i));
    w_41(i)=((a*e)-(b*d))/((c*e)-(b*f));
    w_31(i)=((a*f)-(c*d))/((b*f)-(c*e));
    end
    
    %Cálculo de las veclocidades y la aceleración de B
 
    vB(i)=w_21*l_2; %como w_21 es constante y su distanica al CIR también, 
                    %su velocidad es constante
                    
    vC(i)=w_41(i)*l_4;%usando la fórmula Vc=w_41·R
    
    aB(i)=(w_21^2)*l_2;%siendo la vB cte, solo tendrá componente normal
                       %siendo su fórmula aB=w_21^2·R
end
%% Cálculo de las aceleraciones angulares
%Apoyándome en el planteamiento de las ecuaiones que adjunté
for i=1:1:N
    if  phi_2(i)<alpha_1 && phi_2(i)>=(pi()/2)
   
    a=(w_21^2)*l_2*cos(phi_2(i))-((w_31(i))^2)*l_3*cos(phi_3(i))-((w_41(i))^2)*l_4*cos(phi_4(i));
    b=l_3*sin(phi_3(i));
    c=-l_4*sin(phi_4(i));
    d=-(w_21^2)*l_2*sin(phi_2(i))+((w_31(i))^2)*l_3*sin(phi_3(i))+((w_41(i))^2)*l_4*sin(phi_4(i));
    e=l_3*cos(phi_3(i));
    f=-l_4*cos(phi_4(i));
    dw_41dt(i)=((a*e)-(b*d))/((c*e)-(b*f));
    dw_31dt(i)=((a*f)-(c*d))/((b*f)-(c*e));
    A=-dw_41dt(i)*l_4*sin(phi_4(i))+(w_41(i)^2)*l_4*cos(phi_4(i));
    B=-dw_41dt(i)*l_4*cos(phi_4(i))-(w_41(i)^2)*l_4*sin(phi_4(i));
    aC(i)=sqrt((A^2)+(B^2));
    
    elseif phi_2(i)<(pi()/2) && phi_2(i)>=alpha_2
        
    a=-(w_21^2)*l_2*cos(phi_2(i))-((w_31(i))^2)*l_3*cos(phi_3(i))-((w_41(i))^2)*l_4*cos(phi_4(i));
    b=l_3*sin(phi_3(i));
    c=l_4*sin(phi_4(i));
    d=-(w_21^2)*l_2*sin(phi_2(i))-((w_31(i))^2)*l_3*sin(phi_3(i))+((w_41(i))^2)*l_4*sin(phi_4(i));
    e=-l_3*cos(phi_3(i));
    f=l_4*cos(phi_4(i));
    dw_41dt(i)=((a*e)-(b*d))/((c*e)-(b*f));
    dw_31dt(i)=((a*f)-(c*d))/((b*f)-(c*e));
    A=dw_41dt(i)*l_4*sin(phi_4(i))+(w_41(i)^2)*l_4*cos(phi_4(i));
    B=dw_41dt(i)*l_4*cos(phi_4(i))-(w_41(i)^2)*l_4*sin(phi_4(i));
    aC(i)=sqrt((A^2)+(B^2));
    
    elseif phi_2(i)<alpha_2 && phi_2(i)>=0

    a=-(w_21^2)*l_2*cos(phi_2(i))-((w_31(i))^2)*l_3*cos(phi_3(i))-((w_41(i))^2)*l_4*cos(phi_4(i));
    b=l_3*sin(phi_3(i));
    c=-l_4*sin(phi_4(i));
    d=-(w_21^2)*l_2*sin(phi_2(i))-((w_31(i))^2)*l_3*sin(phi_3(i))+((w_41(i))^2)*l_4*sin(phi_4(i));
    e=-l_3*cos(phi_3(i));
    f=-l_4*cos(phi_4(i));
    dw_41dt(i)=((a*e)-(b*d))/((c*e)-(b*f));
    dw_31dt(i)=((a*f)-(c*d))/((b*f)-(c*e));
    A=-dw_41dt(i)*l_4*sin(phi_4(i))+(w_41(i)^2)*l_4*cos(phi_4(i));
    B=-dw_41dt(i)*l_4*cos(phi_4(i))-(w_41(i)^2)*l_4*sin(phi_4(i));
    aC(i)=sqrt((A^2)+(B^2));
    
    end
end
%% Cálculo de los ángulos de la barras 2,3 y 4 en sentido antihorario y respecto al eje X

%Paso de radianes a sexagesimales
phi_2=(180/pi())*phi_2;
phi_3=(180/pi())*phi_3;
phi_4=(180/pi())*phi_4;
alpha_1=(180/pi())*alpha_1;
alpha_2=(180/pi())*alpha_2;

for i=1:1:N
    if  phi_2(i)<alpha_1 && phi_2(i)>=90
    phi_3(i)=abs(phi_3(i)-180);
    phi_4(i)=abs(phi_4(i)-180);
    elseif phi_2(i)<90 && phi_2(i)>=alpha_2
    phi_4(i)=abs(phi_4(i)-180);
    elseif phi_2(i)<alpha_2 && phi_2(i)>=0    
    phi_4(i)=abs(phi_4(i)-180);
    end
end
%% Añadición del instante inicial

    %Creación de nuevos vectores de dimensión N+1 conteniendo ahora la
    %condición inicial
    It=zeros(1,N+1);
    Iphi_3=zeros(1,N+1);
    Iphi_4=zeros(1,N+1);
    Iw_31=zeros(1,N+1);
    Iw_41=zeros(1,N+1);
    Idw_31dt=zeros(1,N+1);
    Idw_41dt=zeros(1,N+1);
    IvB=zeros(1,N+1);
    IvC=zeros(1,N+1);
    IaB=zeros(1,N+1);
    IaC=zeros(1,N+1);
    
    %Sabiendo que w_21 en el instante inicial es 1 (rad/s)
    %Condición inicial (mismos valores cuando phi_2 = 90 ya calculados) 
    Iphi_3(1)=180;
    Iphi_4(1)=135;
    Iw_31(1)=1;
    Iw_41(1)=1;
    Idw_31dt(1)=2;
    Idw_41dt(1)=2;
    IvB(1)=0.1;
    IvC(1)=0.141;
    IaB(1)=0.1;
    IaC(1)=0.316;
    
for i=1:1:N
    
    It(i+1)=t(i);
    Iphi_3(i+1)=phi_3(i);
    Iphi_4(i+1)=phi_4(i);
    Iw_31(i+1)=w_31(i);
    Iw_41(i+1)=w_41(i);
    Idw_31dt(i+1)=dw_31dt(i);
    Idw_41dt(i+1)=dw_41dt(i);
    IvB(i+1)=vB(i);
    IvC(i+1)=vC(i);
    IaB(i+1)= aB(i);
    IaC(i+1)=aC(i);
end

%% Gráfica de cada variable en función del tiempo

subplot(5,2,1)
plot(It,Iphi_3,'b','linewidth',1.3)
xlabel('tiempo (s)')
ylabel('grados sexagesimales')
legend({'phi_3'},'Location','northeast')
grid on; grid minor;
title('Ángulos de la barra 3 con respecto al eje X en función de t ')
ylim([0.5,200])

subplot(5,2,2)
plot(It,Iphi_4,'r','linewidth',1.3)
xlabel('tiempo (s)')
ylabel('grados sexagesimales')
legend({'phi_4'},'Location','northeast')
grid on; grid minor;
title('Ángulos de la barra 4 con respecto al eje X en función de t ')
ylim([100,155])

subplot(5,2,3)
plot(It,Iw_31,'b','linewidth',1.3)
xlabel('tiempo (s)')
ylabel('velocidad angular (rad/s)')
legend({'w_3'},'Location','northeast')
grid on; grid minor;
title('Velocidad de la barra 3 en función de t respecto al sistema fijo')

subplot(5,2,4)
plot(It,Iw_41,'r','linewidth',1.3)
xlabel('tiempo (s)')
ylabel('velocidad angular (rad/s)')
legend({'w_4'},'Location','northeast')
grid on; grid minor;
title('Velocidad de la barra 4 en función de t respecto al sistema fijo')

subplot(5,2,5)
plot(It,Idw_31dt,'b','linewidth',1.3)
xlabel('tiempo (s)')
ylabel('aceleración angular (rad/s^2)')
legend({'d(w_3)/dt'},'Location','northeast')
grid on; grid minor;
title('Aceleración angular de la barra 3 en función de t respecto al sistema fijo')

subplot(5,2,6)
plot(It,Idw_41dt,'r','linewidth',1.3)
xlabel('tiempo (s)')
ylabel('aceleración angular (rad/s^2)')
legend({'d(w_4)/dt'},'Location','northeast')
grid on; grid minor;
title('Aceleración angular de la barra 3 en función de t respecto al sistema fijo')

subplot(5,2,7)
plot(It,IvB,'b','linewidth',1.3)
xlabel('tiempo (s)')
ylabel('velocidad lineal (m/s)')
legend({'v_B'},'Location','southeast')
grid on; grid minor;
title('Velocidad lineal del punto B en función de t respecto al sistema fijo')
ylim([0,0.15])

subplot(5,2,8)
plot(It,IvC,'r','linewidth',1.3)
xlabel('tiempo (s)')
ylabel('velocidad lineal (m/s)')
legend({'v_c'},'Location','northeast')
grid on; grid minor;
title('Velocidad lineal del punto C en función de t respecto al sistema fijo')

subplot(5,2,9)
plot(It,IaB,'b','linewidth',1.3)
xlabel('tiempo (s)')
ylabel('velocidad lineal (m/s^2)')
legend({'a_B'},'Location','southeast')
grid on; grid minor;
title('Aceleración del punto B en función de t respecto al sistema fijo')
ylim([0,0.15])

subplot(5,2,10)
plot(It,IaC,'r','linewidth',1.3)
xlabel('tiempo (s)')
ylabel('velocidad lineal (m/s^2)')
legend({'a_c'},'Location','northeast')
grid on; grid minor;
title('Aceleración del punto C en función de t respecto al sistema fijo')