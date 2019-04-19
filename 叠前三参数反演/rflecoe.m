function [coe] = rflecoe(vp,vs,den,ang,mod)

% Note：角度必须是弧度制

m=length(vp);
n=length(ang);
coe=zeros(m,n);

if (strcmp(mod,'zoeppritz'))
    for i=1:n
        alph=ang(i);
        for j=1:m-1
            vp1=vp(j);vp2=vp(j+1);
            vs1=vs(j);vs2=vs(j+1);
            rho1=den(j);rho2=den(j+1);
            
            snell_ratio=sin(alph)/vp1;     % 射线参数，斯奈尔比
            beta=asin(vs1*snell_ratio);
            alph_=asin(vp2*snell_ratio);
            beta_=asin(vs2*snell_ratio);

            %计算A矩阵
            A(1,1)=sin(alph);
            A(1,2)=cos(beta);
            A(1,3)=-1*sin(alph_);
            A(1,4)=cos(beta_);

            A(2,1)=cos(alph);
            A(2,2)=-1*sin(beta);
            A(2,3)=cos(alph_);
            A(2,4)=sin(beta_);

            A(3,1)=cos(2*beta);
            A(3,2)=-1*vs1/vp1*sin(2*beta);
            A(3,3)=-1*(rho2*vp2)/(rho1*vp1)*cos(2*beta_);
            A(3,4)=-1*(rho2*vs2)/(rho1*vp1)*sin(2*beta_);

            A(4,1)=(vs1*vs1)/vp1*sin(2*alph);
            A(4,2)=vs1*cos(2*beta);
            A(4,3)=(rho2*vs2*vs2)/(rho1*vp2)*sin(2*alph_);
            A(4,4)=-1*(rho2*vs2)/rho1*cos(2*beta_);

            %计算B矩阵
            B(1,1)=-1*sin(alph); B(2,1)=cos(alph);  
            B(3,1)=-1*cos(2*beta); 
            B(4,1)=(vs1*vs1)/vp1*sin(2*alph);

            %计算R矩阵
            R=pinv(A)*B;
            coe(j,i)=R(1,1);
        end
    end
end

if (strcmp(mod,'aki_richard'))
    for i=1:n
        alph=ang(i);
        for j=1:m-1
            vp1=vp(j);vp2=vp(j+1);
            vs1=vs(j);vs2=vs(j+1);
            rho1=den(j);rho2=den(j+1);
            alph_=asin(vp2/vp1*sin(alph));
            alph=0.5*(alph+alph_);
            dvp=vp2-vp1;mvp=0.5*(vp2+vp1);
            dvs=vs2-vs1;mvs=0.5*(vs2+vs1);
            drho=rho2-rho1;mrho=0.5*(rho2+rho1);
            a_alph=0.5*(1+tan(alph)*tan(alph));
            a_beta=-4*mvs*mvs*sin(alph)*sin(alph)/(mvp*mvp);
            a_rho=0.5*(1+a_beta);
            
            coe(j,i)=a_alph*dvp/mvp+a_beta*dvs/mvs+a_rho*drho/mrho;
        end
    end
end

if (strcmp(mod,'shuey'))
    for i=1:n
        alph=ang(i);
        for j=1:m-1
            vp1=vp(j);vp2=vp(j+1);
            vs1=vs(j);vs2=vs(j+1);
            rho1=den(j);rho2=den(j+1);
            alph_=asin(vp2/vp1*sin(alph));
            alph=0.5*(alph+alph_);
            dvp=vp2-vp1;mvp=0.5*(vp2+vp1);
            dvs=vs2-vs1;mvs=0.5*(vs2+vs1);
            drho=rho2-rho1;mrho=0.5*(rho2+rho1);
            P=0.5*(dvp/mvp+drho/mrho);
            G=0.5*dvp/mvp-2*((mvs/mvp)^2)*(drho/mrho)-4*((mvs/mvp)^2)*(dvs/mvs);
            C=0.5*dvp/mvp;
            
%             coe(j,i)=P+G*sin(alph)*sin(alph);          %shuey 二项式
            coe(j,i)=P+G*sin(alph)*sin(alph)+C*(tan(alph)*tan(alph)-sin(alph)*sin(alph));  %shuey 三项式
%             a_alph=0.5*(1+tan(alph)*tan(alph));
%             a_beta=-4*mvs*mvs*sin(alph)*sin(alph)/(mvp*mvp);
%             a_rho=0.5*(1+a_beta);
%             
%             coe(j,i)=a_alph*dvp/mvp+a_beta*dvs/mvs+a_rho*drho/mrho;
        end
    end
end






           
