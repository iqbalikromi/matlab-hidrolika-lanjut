% UTS HIDRAULIKA II
% AUTHOR: AKHMAD IQBAL IKROMI
% COPY-PASTE ISN'T ALLOWED, READ AND LEARN ONLY

clear
clf

yn=3.917; %menggunakan data no 1
b=15;

L=10000; %m
dt=5; %detik
dx=500;
g=9.81;
m=0;
n=0.01; %frictionless
So=0.0002;
elv_hulu=1; %asumsi elevasi di hulu

periode_pasut= 4000; %detik, simulasi minimal 3kali periode

time=3*periode_pasut;
tm=time/dt;
im=L/dx+1;
out=1300;
tpl=time/out;

%initial condition
%polutan
for i=1:im
     
    Q(i)=120;
    A(i)=b*yn;
    y(i)=yn;
    x(i)=(i-1)*dx;
    z(i)=elv_hulu-x(i)*So;

% UTS HIDRAULIKA II
% AUTHOR: AKHMAD IQBAL IKROMI
% COPY-PASTE ISN'T ALLOWED, READ AND LEARN ONLY

end
% plot(x,z,'color','black','LineWidth',3,'DisplayName','Bed Elevation')
% hold on

%program
for t=1:tm
    ts(t)=t*dt;

    %kontinuitas
    for i=2:im-1
        Ab(i)=A(i)-dt/2/dx*(Q(i+1)-Q(i-1));
        yb(i)=Ab(i)/b;
    end
    %momentum
    for i=2:im-1
        u=Q(i)/A(i);
        y(i)=A(i)/b;
        R=(b*y(i))/(b+2*y(i));
        Sf=n^2*abs(u)*u/R^(4/3);
        bro1=dt/2/dx*(Q(i+1)^2/A(i+1)-Q(i-1)^2/A(i-1));
        bro2=g*A(i)*dt*(So-Sf);
        bro3=g*A(i)*dt/2/dx*(y(i+1)-y(i-1));
        Qb(i)=Q(i)-bro1+bro2-bro3;
    end

% UTS HIDRAULIKA II
% AUTHOR: AKHMAD IQBAL IKROMI
% COPY-PASTE ISN'T ALLOWED, READ AND LEARN ONLY

    %BC hulu
    yb(1)=yb(2);
    Ab(1)=Ab(2);
    Qb(1)=120;
    
    %BC hilir
    Qb(im)=Qb(im-1);
    yb(im)=yn+sin(2*pi*ts(t)/4000);
    Ab(im)=b*yb(im);
   
    A=Ab;
    Q=Qb;
    y=yb;

    %FILTER HANSEN untuk hidrodinamik
    Afil=A;
    Qfil=Q;
    Cofil=0.99;
    for i=2:im-1
        A(i)=Cofil*Afil(i)+(1-Cofil)*(Afil(i-1)+Afil(i+1))/2;
        Q(i)=Cofil*Qfil(i)+(1-Cofil)*(Qfil(i-1)+Qfil(i+1))/2;
    end
       
    for i=8000/dx+1
        yp(t)=y(i);
    end
    for i=6000/dx+1
        yq(t)=y(i);
    end
     for i=4000/dx+1
        yr(t)=y(i);
     end
%     if mod(t*dt,12000)<=0
%         hold on
%         plot(ts,yp,'DisplayName','Pasang Surut di 2km dari hilir'); 
%         plot(ts,yq,'DisplayName','Pasang Surut di 4km dari hilir');   
%         plot(ts,yr,'DisplayName','Pasang Surut di 6km dari hilir');
%         legend(Location="bestoutside")
%         ylabel('Tinggi Muka Air (m)');
%         xlabel('Waktu (s)');
%         title('Hasil Simulasi Pasang Surut','FontSize',16)
%     end
%         if mod(t*dt,720)<=0
% %         hold on
%         yyaxis left
%         txt1=['Muka Air t= ',num2str(t*dt/3600), ' jam'];
%         plot(x,(y+z),'blue','DisplayName',txt1);
%         ylim([-2 7])
%         xlim([0 10000])
%         ylabel('Tinggi Muka Air (m)')
%         xlabel('Jarak (m)');
%         yyaxis right
%         ylim([-2 7])
%         title('Hasil Simulasi Pasang Surut','FontSize',16)
%         legend("show" ,Location="bestoutside")
%         pause(0.1);
%         hold off
%       
%     end
% end

% UTS HIDRAULIKA II
% AUTHOR: AKHMAD IQBAL IKROMI
% COPY-PASTE ISN'T ALLOWED, READ AND LEARN ONLY



    %plotting ANIMASI 
    if mod(t*dt,tpl)<=0
        p1=plot(y+z,'blue');
        hold on
        p2=plot(z,'red');
        axis([2000/dx im -2 10]);
        title('Simulasi Pasut');
        legend([p1,p2],'Muka Air','Bed Elv');
        xlabel('i*dx')
        ylabel('muka air');
        pause(0.1);
        hold off;
        end
end





