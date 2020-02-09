clc
clear
% close all
KK=[-1.0801    -0.7715  -0.4629  -0.1543  0.1543 0.4629  0.7715  1.0801];

delta=0;
beta_skew=0;


zel=[-1.0801  -0.7715  -0.4629 -0.1543  0.1543  0.4629  0.7715  1.0801];
sefr=zeros(1,8);
alpha=1.7;
SNR_inx=1;
EQ=[-2.5:0.05:2.5];

SNR=-2:2:8;
kappa=5;
% gamma_g=((inv(10^(SNR/10)))*0.5)*(1/(kappa+1));
% gamma_s=kappa*gamma_g;


% gamma_s=0.0791;
% gamma_g=0.0395;
% 
% noise_total=gamma_s+gamma_g;


loop=1000;

%% ACTUAL LLR
%%  Bit 3

Gauss=[];
Gauss_neg=[];

for v=1:length(SNR)
    kk=0;
    for KJ=[-2.5:0.05:2.5];
        gamma_s=((inv(10^(SNR(v)/10)))*0.5)*(1/(kappa+1));
        gamma_g=kappa*gamma_s;
        kk=kk+1;
        Dist_1=-0.7715;
        Dist_2=-0.4629;
        Dist_3=0.46291;
        Dist_4=0.7715;
        Dist_5=-1.0801;
        Dist_6=-0.1543;
        Dist_7=0.1543;
        Dist_8=1.0801;
        f_1=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_1)*x));
        f_2=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_2)*x));
        f_3=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_3)*x));
        f_4=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_4)*x));
        f_5=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_5)*x));
        f_6=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_6)*x));
        f_7=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_7)*x));
        f_8=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_8)*x));
        Gauss(1,kk) = integral(f_1,0,1000);
        Gauss(2,kk) = integral(f_2,0,1000);
        Gauss(3,kk) = integral(f_3,0,1000);
        Gauss(4,kk) = integral(f_4,0,1000);
        Gauss_neg(1,kk) = integral(f_5,0,1000);
        Gauss_neg(2,kk) = integral(f_6,0,1000);
        Gauss_neg(3,kk) = integral(f_7,0,1000);
        Gauss_neg(4,kk) = integral(f_8,0,1000);
    end
    LLR_Actual_3_mat(v,:)=log((sum(Gauss))./(sum(Gauss_neg)));
end
%%  Bit 2
kk=0;
Gauss=[];
Gauss_neg=[];
for v=1:length(SNR)
    v
    kk=0;
    for KJ=[-2.5:0.05:2.5];
        gamma_s=((inv(10^(SNR(v)/10)))*0.5)*(1/(kappa+1));
        gamma_g=kappa*gamma_s;
        kk=kk+1;
        Dist_1=0.1543 ;
        Dist_2=0.4629;
        Dist_3=-0.1543;
        Dist_4=-0.4629;
        Dist_5=0.7715;
        Dist_6=1.0801;
        Dist_7=-0.7715;
        Dist_8=-1.0801;
        f_1=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_1)*x));
        f_2=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_2)*x));
        f_3=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_3)*x));
        f_4=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_4)*x));
        f_5=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_5)*x));
        f_6=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_6)*x));
        f_7=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_7)*x));
        f_8=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_8)*x));
        Gauss(1,kk) = integral(f_1,0,1000);
        Gauss(2,kk) = integral(f_2,0,1000);
        Gauss(3,kk) = integral(f_3,0,1000);
        Gauss(4,kk) = integral(f_4,0,1000);
        Gauss_neg(1,kk) = integral(f_5,0,1000);
        Gauss_neg(2,kk) = integral(f_6,0,1000);
        Gauss_neg(3,kk) = integral(f_7,0,1000);
        Gauss_neg(4,kk) = integral(f_8,0,1000);
    end
    LLR_Actual_2_mat(v,:)=log((sum(Gauss))./(sum(Gauss_neg)));
end


%%  Bit 1
kk=0;
Gauss=[];
Gauss_neg=[];

for v=1:length(SNR)
    kk=0;
    for KJ=[-2.5:0.05:2.5];
        
        kk=kk+1;
        gamma_s=((inv(10^(SNR(v)/10)))*0.5)*(1/(kappa+1));
        gamma_g=kappa*gamma_s;
        Dist_1=0.1543 ;
        Dist_2=0.4629;
        Dist_3=0.7715;
        Dist_4=1.0801;
        Dist_5=-0.1543 ;
        Dist_6=-0.4629 ;
        Dist_7=-0.7715;
        Dist_8=-1.0801;
        f_1=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_1)*x));
        f_2=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_2)*x));
        f_3=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_3)*x));
        f_4=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_4)*x));
        f_5=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_5)*x));
        f_6=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_6)*x));
        f_7=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_7)*x));
        f_8=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_8)*x));
        Gauss(1,kk) = integral(f_1,0,1000);
        Gauss(2,kk) = integral(f_2,0,1000);
        Gauss(3,kk) = integral(f_3,0,1000);
        Gauss(4,kk) = integral(f_4,0,1000);
        Gauss_neg(1,kk) = integral(f_5,0,1000);
        Gauss_neg(2,kk) = integral(f_6,0,1000);
        Gauss_neg(3,kk) = integral(f_7,0,1000);
        Gauss_neg(4,kk) = integral(f_8,0,1000);
    end
    LLR_Actual_1_mat(v,:)=log((sum(Gauss))./(sum(Gauss_neg)));
end


% alpha=1.5;
SNR=[-2:2:8]
for v=1:length(SNR)
    gamma_s=((inv(10^(SNR(v)/10)))*0.5)*(1/(kappa+1));
    gamma_g=kappa*gamma_s;
    sigma_n=2*gamma_g;
    A_dis=makedist('Stable','alpha',alpha,'beta',beta_skew,'gam',gamma_s,'delta',0);
    sum_e=0;
    sum_prop=0;
    sum_th=0;
    sum_dif=0;
    sum_bit=0;
    for i=1:loop
        data = (randi([0 1],6,1));
        Sym=qammod(data,64, 'InputType','bit', 'UnitAveragePower',true);
        
        Gaus_real=normrnd(0,sqrt(sigma_n));
        Gaus_im=normrnd(0,sqrt(sigma_n));
        Alpha_real=random(A_dis);
        Alpha_img=random(A_dis);
        
        Recieved_real=real(Sym)+Alpha_real+Gaus_real;
        Recieved_im=imag(Sym)+Alpha_img+Gaus_im;
        rec=Recieved_real+ j*Recieved_im;
        LLR_qam=qamdemod(rec, 64,'OutputType','llr', 'NoiseVariance',sigma_n, 'UnitAveragePower',true);
        QAM_ROW=zeros(1,6);
        QAM_ROW((LLR_qam>=0))=1;
        MAT_1(i,:)=QAM_ROW;
        bit_qam=qamdemod(rec, 64,'OutputType','bit','UnitAveragePower',true);
        DEQAM_out(i,:)=LLR_qam;
        %     bit_qam=qamdemod(rec, 64,'OutputType','llr');
        Bit_m=zeros(1,6);
        Bit_m((LLR_qam<0))=1;
        [E1,E2]=biterr(Bit_m',data);
        [G1,G2]=biterr(bit_qam,data);
        sum_e=sum_e+E1;  
        sum_bit=sum_bit+G1;
        LLR_real_1(i)=LLR_1(Recieved_real,EQ,gamma_g,gamma_s,alpha,v);
        LLR_real_2(i)=LLR_2(Recieved_real,EQ,gamma_g,gamma_s,alpha,v);
        LLR_real_3(i)=LLR_3(Recieved_real,EQ,gamma_g,gamma_s,alpha,v);
        LLR_im_1(i)=LLR_1(Recieved_im,EQ,gamma_g,gamma_s,alpha,v);
        LLR_im_2(i)=LLR_2(Recieved_im,EQ,gamma_g,gamma_s,alpha,v);
        LLR_im_3(i)=LLR_3(Recieved_im,EQ,gamma_g,gamma_s,alpha,v);
        LLR_est_t=[ -LLR_real_1(i)  -LLR_real_2(i)   -LLR_real_3(i)   LLR_im_1(i) -LLR_im_2(i)  -LLR_im_3(i)];
        Bit_prop=zeros(1,6);
        Bit_prop((LLR_est_t<0))=1;
        [W1,W2]=biterr(Bit_prop',data);
        sum_prop=sum_prop+W1;
        
        
        
        
        [e1,e2]=min(abs(EQ-Recieved_real));
        LLR_Actual_1=LLR_Actual_1_mat(v,:);
        LLR_Actual_2=LLR_Actual_2_mat(v,:);
        LLR_Actual_3=LLR_Actual_3_mat(v,:);
        ACT_real_1(i)=LLR_Actual_1(e2);
        ACT_real_2(i)=LLR_Actual_2(e2);
        ACT_real_3(i)=LLR_Actual_3(e2);
        
        [c1,c2]=min(abs(EQ-Recieved_im));
        ACT_img_1(i)=LLR_Actual_1(c2);
        ACT_img_2(i)=LLR_Actual_2(c2);
        ACT_img_3(i)=LLR_Actual_3(c2);
        
        LLR_the_t=[ -ACT_real_1(i)  -ACT_real_2(i)   -ACT_real_3(i)   ACT_img_1(i) -ACT_img_2(i)  -ACT_img_3(i)];
        Bit_th=zeros(1,6);
        Bit_th((LLR_the_t>=0))=1;
        MAT_2(i,:)=Bit_th;
        B_TH=zeros(1,6);
        B_TH((LLR_the_t<0))=1;
        [Z1,Z2]=biterr(B_TH',data);
        sum_th=sum_th+Z1;
        
        [TE_1,TE_2]=biterr(B_TH',Bit_prop');
        sum_dif=sum_dif+TE_1;
        x_i(i)=Recieved_real;
        y_i(i)=Recieved_im;
    end
    [J1,J2]=sort(x_i);
    [J3,J4]=sort(y_i);
%     figure
%     plot(J3,LLR_im_1(J4),'r--');hold on;plot(J3,ACT_img_1(J4),'b--');hold on;plot(J1, DEQAM_out(:,4),'g--')
%     hold on; plot(J3,LLR_im_2(J4),'r-v');hold on;plot(J3,ACT_img_2(J4),'b-v');hold on;plot(J1, DEQAM_out(:,5),'g--')
%     hold on; plot(J3,LLR_im_3(J4),'r-o');hold on;plot(J3,ACT_img_3(J4),'b-o');hold on;plot(J1, DEQAM_out(:,6),'g--')
%     
%     figure
%     plot(J1,LLR_real_1(J2),'r--');hold on;plot(J1,ACT_real_1(J2),'b--');hold on;plot(J1, DEQAM_out(:,1),'g--')
%     hold on ; plot(J1,LLR_real_2(J2),'r-v');hold on;plot(J1,ACT_real_2(J2),'b-v');hold on;plot(J1, DEQAM_out(:,2),'g-v')
%     hold on; plot(J1,1*LLR_real_3(J2),'r-o');hold on;plot(J1,1*ACT_real_3(J2),'b-o');hold on;plot(J1, DEQAM_out(:,3),'g-o')
    err_rate(v)=sum_e/(loop*6);
    err_rate_prop(v)=sum_prop/(loop*6);
    err_rate_th(v)=sum_th/(loop*6);
    err_rate_diff(v)=sum_dif/(loop*6);
    err_bit_qam(v)=sum_bit/(loop*6);
    
end

figure
plot(SNR,err_rate,'r');hold on;plot(SNR,err_rate_prop,'g');hold on;plot(SNR,err_rate_th,'b');hold on;hold on;plot(SNR,err_rate_diff,'k');hold on;
plot(SNR, err_bit_qam,'m')
