clear
clf
clc

time=100;
dt=0.06;
loop=time/dt;
output=time/10;
L=700; %m
b=400; %m

dx=10;
dy=10;


im=L/dx+1;
jm=b/dy+1;

Pe=2;
%initial condition
for i=1:im
    for j=1:jm
        x(i)=(i-1)*dx;
        y(j)=(j-1)*dy;
        u(i,j)=0.5;
        v(i,j)=0;
        h(i,j)=1;
        c(i,j)=0;

    end
end

ispil=100/dx+1;
jspil=200/dx+1;
c(ispil,jspil)=1;


%engine
for t=1:loop
    for i=1:im-2
        for j=1:jm-2
            ac=u(i,j)*dx/Pe;
            %predictor

            suku1p=(u(i+1,j)*h(i+1,j)*c(i+1,j)-u(i,j)*h(i,j)*c(i,j))/dx;
            suku2p=(v(i,j+1)*h(i,j+1)*c(i,j+1)-v(i,j)*h(i,j)*c(i,j))/dy;
            suku3p=ac*h(i,j)*(c(i+2,j)-2*c(i+1,j)+c(i,j))/(2*dx)^2;
            suku4p=ac*h(i,j)*(c(i,j+2)-2*c(i,j+1)+c(i,j))/(2*dy)^2;
            cp(i,j)=(h(i,j)*c(i,j)+dt*(suku3p+suku4p-suku1p-suku2p))/h(i,j);
        end
    end
    for i=3:im-2
        for j=3:jm-2
            %corrector-final
            chalf=0.5*(c(i,j)+cp(i,j));
            suku1c=(u(i,j)*h(i,j)*c(i,j)-u(i-1,j)*h(i-1,j)*c(i-1,j))/dx;
            suku2c=(v(i,j)*h(i,j)*c(i,j)-v(i,j-1)*h(i,j-1)*c(i,j-1))/dy;
            suku3c=ac*h(i,j)*(c(i,j)-2*c(i-1,j)+c(i-2,j))/(2*dx)^2;
            suku4c=ac*h(i,j)*(c(i,j)-2*c(i,j-1)+c(i,j-2))/(2*dy)^2;
            cn(i,j)=(h(i,j)*chalf+dt*(suku3c+suku4c-suku1c-suku2c))/h(i,j);
        end
    end

    %syarat batas
            %bagian ujung kiri ketika j=3
    for j=3:jm-2
        cn(2,j)=cn(3,j);
        %bagian ujung kanan ketika i=im-1
        cn(im-1,j)=cn(im-2,j);
    end
            %dinding kanan i=3
    for i=3:im-2
        cn(i,2)=cn(i,3);
        %dinding kiri
        cn(i,jm-1)=cn(i,jm-2);
    end
        %bagian ujung kiri ketika i=2
    for j=2:jm-1
        cn(1,j)=cn(2,j);
        %bagian ujung kanan ketika i=im
        cn(im,j)=cn(im-1,j);
    end

        %dinding kanan ketiks j=2
    for i=2:im-1
        cn(i,1)=cn(i,2);
        %dinding kiri
        cn(i,jm)=cn(i,jm-1);
    end

    %bagian sudut
    cn(1,1)=(cn(1,2)+cn(2,1))/2;
    cn(1,jm)=(cn(2,jm)+cn(1,jm-1))/2;
    cn(im,1)=(cn(im-1,1)+cn(im,2))/2;
    cn(im,jm)=(cn(im,jm-1)+cn(im-1,jm))/2;

    %filter
    cfl=cn;
    cof=0.98;
    for i=2:im-1
        for j=2:jm-1
            cn(i,j)=cof*cfl(i,j)+(1-cof)*(cfl(i+1,j)+cfl(i,j+1)+cfl(i,j)+cfl(i-1,j)+cfl(i,j-1))/5;
        end
    end

    %loop
    c=cn;
    if time==10
        zz=0.25;
        zl=0.15;

    elseif time==100
        zz=0.025;
        zl=0.015;

    else
        zz=0.005;
        zl=0.0026;
    end
 
    
    %plotting
    if mod(t*dt,output)==0
        mesh(y,x,c);
        shading("interp");
        axis equal
        set(gca,'DataAspectRatio',[350 400 zz])
%         axis vis3d
        zlim([0 zl]);
%         xlim([0 y(jm)]);
%         ylim([0 x(im)]);
        view(45,30);
        colormap('turbo')
        colorbar
        pause(1);
        
      
    end

end

