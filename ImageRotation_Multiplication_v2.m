% Author: Pushpendra singh
% Paper ref.: On the hypercomplex numbers and normed division algebra of all dimensions: A unified multiplication

% Generate Fig. 1, 8, 9, and 10 of the paper (simulations)

close all;
clear all;
format long;
clc;

set(0,'defaulttextinterpreter','latex');

if 0 % set 0 to 1 for Fig. 1, OKay
    t=-3*pi/2:0.001:3*pi/2;
    
    sinpi=sp_sinpi(t);
    cospi=sp_cospi(t);
        
    subplot(2,1,1)
    plot(t,cospi,'LineWidth',2)
    ylabel('$\cos\pi(\theta)$')
    grid on
    set(gca,'XTick',-3*pi/2:pi/2:3*pi/2) 
    set(gca,'XTickLabel',{'-3\pi/2','-\pi','-\pi/2', '0', '\pi/2','\pi', '3\pi/2'})

    subplot(2,1,2)
    plot(t,sinpi,'LineWidth',2)
    ylabel('$\sin\pi(\theta)$')
    grid on

    set(gca,'XTick',-3*pi/2:pi/2:3*pi/2) 
    set(gca,'XTickLabel',{'-3\pi/2','-\pi','-\pi/2', '0', '\pi/2','\pi', '3\pi/2'})

    %tmpOrthogonal=sum(sinpi.*cospi);    
    tt=1;    
end



if 0 % set 0 to 1 for Fig. 8, Okay

    load("xyzPoints");
    ptCloud = pointCloud(xyzPoints);
    subplot(2,2,1)
    pcshow(ptCloud)
    
    NewIm1=PointCloudImageMultiplication(xyzPoints, xyzPoints); 
    
    ptCloud2 = pointCloud(NewIm1);
    subplot(2,2,2)
    pcshow(ptCloud2) 
    
    
    [X,Y,Z] = sphere(100);
    loc1 = [X(:),Y(:),Z(:)];
    loc2 = loc1/2;
    %loc3 = loc1/2;
    ptCloud = pointCloud([loc1;loc2]);
    subplot(2,2,3)
    pcshow(ptCloud)
    %title('Point Cloud')
    xyzPoints1=ptCloud.Location;
    
    NewIm1=PointCloudImageMultiplication(xyzPoints1, xyzPoints1); 
    
    ptCloud2 = pointCloud(NewIm1);
    subplot(2,2,4)
    pcshow(ptCloud2)
    
    set(gcf,'color','w');
    
    tt=1;    
end



if 0 %set 0 to 1 for Fig. 9, Okay
    load("xyzPoints");
    ptCloud = pointCloud(xyzPoints);
    
    figure
    phi=linspace(0,pi,25);
    tmp=xyzPoints;
    for i=1:length(phi)
        rot=[1,0,phi(i)];
        NewIm=PointCloudRotation(xyzPoints, rot);  
        ptCloud1 = pointCloud(NewIm);
        subplot(5,5,i)
        pcshow(ptCloud1)
        %title('(i)')
        i
    end
    set(gcf,'color','w');
    tt=1; 
end


if 0 % set 0 to 1 for 3D Fig. 10 (left), cube OKay
    array_width=51; %must be odd
    [X,Y,Z,outer_cube]=create_cube(array_width);
    figure
    loc1=([X(:),Y(:),Z(:)]);
    %trisurf(k,X(:),Y(:),Z(:),outer_cube,'EdgeColor','none')
    ptCloud = pointCloud([loc1]);
    %subplot(2,2,3)
    pcshow(ptCloud)
    
    tmp=loc1; 
    
    for i=1:25
        if (i==1)
            rot=[1,0,0];
        else        
            rot=[1,2*pi/24,0]; % Fig 10 (left)
            %rot=[1,0,pi/24]; % Fig 10 (Right)
        end
    NewIm=PointCloudRotation(tmp, rot);  
    tmp=NewIm;
    ptCloud1 = pointCloud(NewIm);
    subplot(5,5,i)
    pcshow(ptCloud1)
    %title('(i)')
    i
    end
    set(gcf,'color','w');
    
    tt=1;  
end

if 1 % set 0 to 1 for 3D Fig. 10 (right), cube OKay
    array_width=51; %must be odd
    [X,Y,Z,outer_cube]=create_cube(array_width);
    figure
    loc1=([X(:),Y(:),Z(:)]);
    %trisurf(k,X(:),Y(:),Z(:),outer_cube,'EdgeColor','none')
    ptCloud = pointCloud([loc1]);
    %subplot(2,2,3)
    pcshow(ptCloud)
    
    tmp=loc1; 
    
    for i=1:25
        if (i==1)
            rot=[1,0,0];
        else        
            %rot=[1,2*pi/24,0]; % Fig 10 (left)
            rot=[1,0,pi/24]; % Fig 10 (Right)
        end
    NewIm=PointCloudRotation(tmp, rot);  
    tmp=NewIm;
    ptCloud1 = pointCloud(NewIm);
    subplot(5,5,i)
    pcshow(ptCloud1)
    %title('(i)')
    i
    end
    set(gcf,'color','w');
    
    tt=1;  
end

function [NewIm]=PointCloudImageMultiplication(g1,g2)
    [R,C]=size(g1);
    NewIm=zeros(R,C);
    for i=1:R    
        NewIm(i,:)=multiplicationMD(g1(i,:),g2(i,:));    
    end
end

function [mult_CCS]=multiplicationMD(g1,g2) % pixel by pixel
    %r=g1(1)*g2(1);        
    %theta=mod(g1(2)+g2(2), 2*pi); % b = mod(a,m), 2pi modulo
    %phi=mod((g1(3)+g2(3)), 2*pi);  % 2pi modulo phi [0, 2pi) 
    %sp=[r;theta;phi];

    N=length(g1);% dimentions
    %N=3;
    
    g1=g1+10^(-15); % ensure that a!=0
    g2=g2+10^(-15); % ensure that b!=0

    g1=Cart2Shperical_phi(g1);
    g2=Cart2Shperical_phi(g2);

    
    mult_SCS=zeros(N,1);
    r=g1(1)*g2(1);
    
    for i=2:2
        theta(i-1)=mod(g1(i)+g2(i), 2*pi); % b = mod(a,m), 2pi modulo         
    end
    for i=3:N
        tmp=g1(i)+g2(i); % ensure phi\in(-pi/2, pi/2)
        if(tmp>pi/2)
           tmp=tmp-pi; 
        end
        if(tmp<-pi/2)
           tmp=tmp+pi; 
        end
        theta(i-1)=tmp;       
    end
    mult_SCS(1,1)=r;
    mult_SCS(2:end,1)=theta;  

    mult_CCS=Shperical2Cart_phi(mult_SCS);
    %mult_CCS=Shperical2Cart_3D(mult_SCS);
    
end



function [NewIm]=ImageMultiplication(g1,g2)
[R,C,D]=size(g1);
NewIm=zeros(R,C,D);
for i=1:R
    for j=1:C
        NewIm(i,j,:)=multiplication(g1(i,j,:),g2(i,j,:));
    end
end
end


function [NewIm]=ImageRotation(g1, rot)
[R,C,D]=size(g1);
NewIm=zeros(R,C,D);

for i=1:R
    for j=1:C
        NewIm(i,j,:)=ImgRotation(g1(i,j,:),rot);
    end
end
end

function [NewIm]=PointCloudRotation(g1, rot)
    [R,C]=size(g1);
    NewIm=zeros(R,C);
    for i=1:R
        NewIm(i,:)=ImgRotation(g1(i,:),rot);        
    end
    tt=1;
end

function [mult_CCS]=ImgRotation(g1,rot) % pixel by pixel
    %r=g1(1)*g2(1);        
    %theta=mod(g1(2)+g2(2), 2*pi); % b = mod(a,m), 2pi modulo
    %phi=mod((g1(3)+g2(3)), 2*pi);  % 2pi modulo phi [0, 2pi) 
    %sp=[r;theta;phi];

    N=length(g1);% dimentions
    %N=3;
    
    %g1=double(reshape(g1,[N,1]));
    %g2=double(reshape(g2,[N,1]))+10^(-15);

    g1=Cart2Shperical_phi(g1);
    g2(:,1)=rot; % it is already in SCS

    
    mult_SCS=zeros(N,1); theta=zeros(N-1,1);
    r=g1(1)*g2(1);
    
    %for i=2:2
        theta(1)=mod(g1(2)+g2(2), 2*pi); % b = mod(a,m), 2pi modulo         
    %end
    for i=3:N
        tmp=g1(i)+g2(i); % ensure phi\in(-pi/2, pi/2)
        if(tmp>pi/2)
           tmp=tmp-pi; 
        end
        if(tmp<-pi/2)
           tmp=tmp+pi; 
        end
        theta(i-1)=tmp;
        %theta(i-1)=mod(g1(i)+g2(i), 2*pi); % b = mod(a,m), 2pi modulo
        %phi=mod((g1(3)+g2(3)), 2*pi);  % 2pi modulo phi [0, 2pi) 
    end
    mult_SCS(1,1)=r;
    mult_SCS(2:end,1)=theta;  

    mult_CCS=Shperical2Cart_phi(mult_SCS);
    %mult_CCS=Shperical2Cart_3D(mult_SCS);

    %mult_CCS=uint8(mult_CCS);
    
    %mult_CCS=reshape(g1,[1,1,3]);

end




function [mult_CCS]=multiplication(g1,g2) % pixel by pixel
    %r=g1(1)*g2(1);        
    %theta=mod(g1(2)+g2(2), 2*pi); % b = mod(a,m), 2pi modulo
    %phi=mod((g1(3)+g2(3)), 2*pi);  % 2pi modulo phi [0, 2pi) 
    %sp=[r;theta;phi];

    N=length(g1);% dimentions
    %N=3;
    
    g1=double(reshape(g1,[N,1]))+10^(-15);
    g2=double(reshape(g2,[N,1]))+10^(-15);

    g1=Cart2Shperical_phi(g1);
    g2=Cart2Shperical_phi(g2);

    
    mult_SCS=zeros(N,1);
    r=g1(1)*g2(1);
    
    for i=2:2
        theta(i-1)=mod(g1(i)+g2(i), 2*pi); % b = mod(a,m), 2pi modulo         
    end
    for i=3:N
        tmp=g1(i)+g2(i); % ensure phi\in(-pi/2, pi/2)
        if(tmp>pi/2)
           tmp=tmp-pi; 
        end
        if(tmp<-pi/2)
           tmp=tmp+pi; 
        end
        theta(i-1)=tmp;
        %theta(i-1)=mod(g1(i)+g2(i), 2*pi); % b = mod(a,m), 2pi modulo
        %phi=mod((g1(3)+g2(3)), 2*pi);  % 2pi modulo phi [0, 2pi) 
    end
    mult_SCS(1,1)=r;
    mult_SCS(2:end,1)=theta;  

    mult_CCS=Shperical2Cart_phi(mult_SCS);
    %mult_CCS=Shperical2Cart_3D(mult_SCS);

    %mult_CCS=uint8(mult_CCS);
    
    %mult_CCS=reshape(g1,[1,1,3]);

end


function [pointin_SCS]=Cart2Shperical_phi(g)
    %a=g(1); b=g(2); c=g(3);
    N=length(g); % Number of dimensions
    pointin_SCS=zeros(N,1); theta=zeros(N-1,1);
    r=sqrt(sum(g.^2));
    pointin_SCS(1,1)=r;
    for i=1:(N-1) 
        if(i==1)
            theta(i)= atan2(g(2),g(1));  % P = atan2(Y,X)             
        else
            %phi=atan2( g(3), sqrt( g(1)*g(1) + g(2)*g(2) ) )\in [-pi/2, pi/2];            
            theta(i)=atan2( g(i+1), sqrt( sum(g(1:i).*g(1:i)) ) ); %            
        end
    end
    pointin_SCS(2:end,1)=theta; 
    %sp=[r;theta;phi]; 
end


function [d]=Shperical2Cart_phi(g)
    N=length(g); d=zeros(N,1);
    r=g(1);
    theta=g(2:end);

    for i=2:N-1 % ensure all phi are in (-pi/2, pi/2)
        phi=theta(i);
        if(phi<-pi/2 || phi>pi/2)
            tmp=mod(phi, pi);
            if(tmp>pi/2)
               tmp=tmp-pi; 
            end
            if(tmp<-pi/2)
               tmp=tmp+pi; 
            end
            theta(i)=tmp;
        end
    end

    tmp=1;
    for i=1:N

        if(i==1)            
            for j=2:(N-1)
                tmp=tmp*cos(theta(j,1));
            end
            tmp=tmp*r*cos(theta(i,1));
            d(1,1)=tmp;
            tmp=1;
        end
        
        if((i>=2) && i<N)            
            for j=i:(N-1)
                tmp=tmp*cos(theta(j,1));
            end
            tmp=tmp*r*sin(theta(i-1,1));
            d(i,1)=tmp;
            tmp=1;
        end
       
        if(i==N)
            tmp=r*sin(theta(i-1,1));  
            d(i,1)=tmp;
        end
               
    end
end

function [cartp]=Shperical2Cart_3D(g)
    r=g(1); theta=g(2); phi=g(3);
    % phi \in (-pi/2,pi/2)
    if(phi<-pi/2 || phi>pi/2)
        tmp=mod(phi, pi);
        if(tmp>pi/2)
           tmp=tmp-pi; 
        end
        if(tmp<-pi/2)
           tmp=tmp+pi; 
        end
        phi=tmp;
    end
    a=r*cos(phi)*cos(theta);
    b=r*cos(phi)*sin(theta);
    c=r*sin(phi);
    cartp=[a;b;c];
end

function [X,Y,Z,cube]=create_cube(array_width)
a=(array_width-1)/2;
a=linspace(-a,a,array_width);
[X,Y,Z]=meshgrid(a);
cube=(X.*Y.*Z);

end

function [x]=sp_sinpi(t)
    t=mod(t,pi);
    k=find(t>=pi/2);
    t(k)=t(k)-pi;
    x=sin(t);
end

function [x]=sp_cospi(t)
    t=mod(t,pi);
    k=find(t>=pi/2);
    t(k)=t(k)-pi;
    x=cos(t);
end

