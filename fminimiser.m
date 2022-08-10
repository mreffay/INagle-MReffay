
function F=fminimiser(X)
%définit le système d'ode du probleme de Laplace
b=X(1);
c=X(2);
global hauteurexp
global largeurexp
global volumeexp
function dY = problemL(s,Y)
dY=zeros(5,1); 
dY(1)=cos(Y(3));
dY(2)=sin(Y(3));
dY(3)=2*b+c*Y(2)-sin(Y(3))/Y(1);
dY(4)=pi*(Y(1))^2*sin(Y(3)); % pourquoi volume c'est ca ? car increment de volume ?
dY(5)=2*pi*Y(1);
%dY1 x
%dY2 z
%dY3 theta
%dY4 Volume
%dY5 Aire
end
%intégration numérique du profil
[s,Y] = ode45(@problemL,[0 0.002],[2;0;0;0;0]); % pourquoi condition initiale c'est 2 ? et pas 0 ?
%condition non mouillante
E=size(Y,1);
for j=1:size(Y,1)
    if Y(j,3)<pi
        E=j;
    end
end
%place l'apex en 0,0
Y(:,1)=Y(:,1)-2;
%calcul valeurs théoriques
hauteurth=Y(E,2);
largeurth=2*max(Y(1:E,1));
volumeth=0;

for k=1:E-1
    volumeth=volumeth+(Y(k+1,2)-Y(k,2))*pi*((Y(k+1,1)+Y(k,1))/2)^2; % pi*x^2*dz
end

%definit la fonction à minimiser 
F=(1-hauteurth/hauteurexp)^2+(1-largeurth/largeurexp)^2+(1-volumeth/volumeexp)^2;
x=Y(:,1);
z=Y(:,2);
z=z-max(z);
z=-z;
[m,I] = min(z);
assignin('base','x',x) 
assignin('base','z',z)
assignin('base','I',I)
% plot(x(1:I),z(1:I),'r')
% plot(-x(1:I),z(1:I),'r')
% Pour tracer juste le profil theorique sans les axes, pour superposition
% plot(x(1:I),z(1:I),'r','LineWidth',5)
% plot(-x(1:I),z(1:I),'r','LineWidth',5)
% set(gca,'visible','off')
% set(gca,'xtick',[])
% set(gca, 'DataAspectRatio', [1 1 1]) 
end 

