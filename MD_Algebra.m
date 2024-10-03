% Author: Pushpendra singh
% Paper ref.: On the hypercomplex numbers and normed division algebras in all dimensions: A unified multiplication

% This code is useful for examining various properties of normed division algebras in N dimensions.

close all;
clear all;
format long;


% Cartesian coordinate system (CCS) and spherical coordinate system (SCS) 

% In general, you can generate N random numbers in the interval 
%(a,b) with the formula r = a + (b-a).*rand(N,1).

a=1; b=-a; N=101; % N is the number of dimentions
NbOfTrial=10000; 
rng('shuffle');
tmpg1=a + (b-a).*rand(N,NbOfTrial); % CCS
rng('shuffle');
tmpg2=a + (b-a).*rand(N,NbOfTrial); % CCS
rng('shuffle');
tmpg3=a + (b-a).*rand(N,NbOfTrial); % CCS


g_bar=zeros(N,1); 
g_minus=zeros(N,1);
r_square_SCS=zeros(N,1);
minus_one_SCS=zeros(N,1); 
minus_one_SCS(1,1)=1; minus_one_SCS(2,1)=pi; minus_one_SCS(3:end,1)=0; 

percent=10;

for i=1:NbOfTrial    
    g1=tmpg1(:,i); % I hypercomplex number of dimention N
    g2=tmpg2(:,i); % II hypercomplex number of dimention N
    g3=tmpg3(:,i); % III hypercomplex number of dimention N
    
    g=g1; % CCS %g_bar=[g1(1); -g1(2); -g1(3)];
    g_bar(1,1)=g(1);
    g_bar(2:end,1)=-g(2:end); %g_minus=[-g1(1); -g1(2); -g1(3)];
    g_minus(1:2)=-g(1:2); g_minus(3:end)=g(3:end); 
    
    
    % we want (-g)g=(-1)g^2; (-g)(-g)=g^2; gg*=r^2; (-g)g*=-r^2
    g_SCS=Cart2Shperical(g);
    g_bar_SCS=Cart2Shperical(g_bar);
    g_minus_SCS=Cart2Shperical(g_minus);
    
    g_zero(:,i)=Shperical2Cart(g_SCS)-g; % ideally zero 
    
    
    r_square_SCS(1)=sum(g.*g); % r^2
    minusr_square_SCS=multiplication_SCS(r_square_SCS,minus_one_SCS); %(-1)r^2
    
    g_square_SCS=multiplication_SCS(g_SCS,g_SCS); % g^2
    minusgxg_SCS=multiplication_SCS(g_minus_SCS,g_SCS); % (-g)g
    
     
    minus_onexgsquare_SCS=multiplication_SCS(minus_one_SCS,g_square_SCS); %(-1)g^2
    
    is_minusgxg_equal2_minusg_square(:,i)=(minusgxg_SCS - minus_onexgsquare_SCS); % (-g)g=(-1)g^2 ?
    
    minusgxminusg_SCS=multiplication_SCS(g_minus_SCS,g_minus_SCS); % (-g)(-g)
    
    is_minusgxminusg_equal2g_square(:,i)=g_square_SCS - minusgxminusg_SCS; % g^2=(-g)(-g) ? should be close to zero
    is_gxbbar_equal2gmod_square(:,i)=multiplication_SCS(g_SCS,g_bar_SCS) - r_square_SCS; % gg*=r^2 should be close to zero
    
    minusgxgbar_SCS=multiplication_SCS(g_minus_SCS,g_bar_SCS); 
    is_minusgxgbar_SCS_equal2minusr_square(:,i)=minusgxgbar_SCS - minusr_square_SCS;% (-g)g*=-r^2 should be close to zero
    
    
    % Convert to SCS   
    g1_SCS=Cart2Shperical(g1);
    g2_SCS=Cart2Shperical(g2);
    g3_SCS=Cart2Shperical(g3);

    g1g2_SCS=multiplication_SCS(g1_SCS,g2_SCS);
    g1g3_SCS=multiplication_SCS(g1_SCS,g3_SCS);
    g2g3_SCS=multiplication_SCS(g2_SCS,g3_SCS);

    % additive associativity g1+(g2+g3)=(g1+g2)+g3; it is obvious and proof is not
    % required

    % mutiplicative associativity g1(g2g3)=(g1g2)g3
    g1mg2g3_SCS=multiplication_SCS(g2g3_SCS,g1_SCS);
    g1g2mg3_SCS=multiplication_SCS(g1g2_SCS,g3_SCS);
    associativeLeftMinusRightSide(:,i)=g1g2mg3_SCS-g1mg2g3_SCS; % ideally should be zero, practically in order of 10^(-15)

    % it is not distributive g1(g2+g3)!=g1g2+g1g3
    %g1g2_CCS=Shperical2Cart(g1g2_SCS);
    %g1g3_CCS=Shperical2Cart(g1g3_SCS);

    %g1xg2plusg3_SCS=multiplication_SCS(g1_SCS, Cart2Shperical(g2+g3));
    %g1xg2plusg3_CCS=Shperical2Cart(g1xg2plusg3_SCS);
    %g1g2pg1g3=g1g2_CCS+g1g3_CCS;

    if(i==NbOfTrial*percent/100)
     percent
     percent=percent+10;
    end

    tt=1;  
end 

    tmp(1)=max(max(abs(associativeLeftMinusRightSide))); % ideally should be zero, practically in order of 10^(-15)
    tmp(2)=max(max(abs(is_minusgxg_equal2_minusg_square))); % (-g)g=-g^2
    tmp(3)=max(max(abs(is_minusgxminusg_equal2g_square)));  % (-g)(-g)=g^2
    tmp(4)=max(max(abs(is_gxbbar_equal2gmod_square)));      % gg*=|g|^2
    %tmp(5)=max(max(abs(is_minusgxgbar_SCS_equal2minusr_square))); % (-g)g*=-r^2
    tmp
    back_forth_CCS2SCS2CCS=max(max(abs(g_zero)))
    
    tt=1;
    
 
function [pointin_SCS]=Cart2Shperical(g)
    %a=g(1); b=g(2); c=g(3);
    N=length(g);
    pointin_SCS=zeros(N,1);
    r=sqrt(sum(g.^2));
    pointin_SCS(1,1)=r;
    theta=zeros(N-1,1);
    for i=1:(N-1) 
        if(i==1)
            theta(i)= atan2(g(2),g(1));  % P = atan2(Y,X)  azimuth angle           
        else
            %atan2( g(3), sqrt( g(1)*g(1) + g(2)*g(2) ) )\in [-pi/2, pi/2] elevation angles;
            theta(i)=atan2( g(i+1), sqrt( sum(g(1:i).*g(1:i)) ) ); 
        end
    end
    pointin_SCS(2:end,1)=theta; 
    %sp=[r;theta;phi]; 
end

function [d]=Shperical2Cart(g)
    N=length(g); d=zeros(N,1);
    r=g(1);
    theta=g(2:end);

    for i=2:N-1
        phi=theta(i);
        if(phi<-pi/2 || phi>pi/2)
            %tmp=mod(phi, pi);
            tmp=phi;
            while(tmp>pi/2)
               tmp=tmp-pi; 
            end
            while(tmp<-pi/2)
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

function [mult_SCS]=multiplication_SCS(g1,g2)
    %r=g1(1)*g2(1);     
    %sp=[r;theta1;theta2];
    N=length(g1);
    mult_SCS=zeros(N,1);
    r=g1(1)*g2(1);
    for i=2:2
        theta(i-1)=mod(g1(i)+g2(i), 2*pi); % b = mod(a,m), 2pi modulo for azimuth angle         
    end
    for i=3:N  % ensure elevation angle \in[-pi/2, pi/2]
        tmp=g1(i)+g2(i);
        while(tmp>pi/2)
           tmp=tmp-pi; 
        end
        while(tmp<-pi/2)
           tmp=tmp+pi; 
        end
        theta(i-1)=tmp;         
    end

    mult_SCS(1,1)=r;
    mult_SCS(2:end,1)=theta;    
end


