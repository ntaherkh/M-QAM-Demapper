clc
clear
close all
KK=[-1.0801    -0.7715  -0.4629  -0.1543  0.1543 0.4629  0.7715  1.0801];

% gamma_s=0.5;
% alpha=1.0;
delta=0;
beta_skew=0;
zel=[-1.0801  -0.7715  -0.4629 -0.1543  0.1543  0.4629  0.7715  1.0801];
sefr=zeros(1,8);
alpha=1.75;
SNR_inx=1;
KJ=[-2.5:0.05:2.5];
% %  KJ(1)=-2.5 ,  KJ(35)=-0.8  KJ(0)=51  KJ(67)=0.8  KJ(101)=2.5 
pos_gauss_start=38;
pos_gauss_end=64;
pos_alpha_left=37;
pos_alpha_right=65;

neg_gauss_start=38;
neg_gauss_end=64;
neg_alpha_left=37;
neg_alpha_right=65;

for SNR=5:5:15
    kappa_inx=1;
    for kappa=1:2:5
        gamma_g=((inv(10^(SNR/10)))*0.5)*(1/(kappa+1));
        gamma_s=kappa*gamma_g;
        
        Dist_1=0.1543 ;
        Dist_2=0.4629;
        Dist_3=0.7715;
        Dist_4=1.0801;
        
        Dist_5=-0.1543 ;
        Dist_6=-0.4629 ;
        Dist_7=-0.7715;
        Dist_8=-1.0801;
        kk=0;
        for KJ=[-2.5:0.05:2.5];
            kk=kk+1;
%             Sf_1=@(x)((exp(-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_1)*x));
%             Sf_2=@(x)((exp(-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_2)*x));
%             Sf_3=@(x)((exp(-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_3)*x));
%             Sf_4=@(x)((exp(-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_4)*x));
%             Sf_5=@(x)((exp(-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_5)*x));
%             Sf_6=@(x)((exp(-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_6)*x));
%             Sf_7=@(x)((exp(-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_7)*x));
%             Sf_8=@(x)((exp(-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_8)*x));
            
            f_1=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_1)*x));
            f_2=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_2)*x));
            f_3=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_3)*x));
            f_4=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_4)*x));
            f_5=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_5)*x));
            f_6=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_6)*x));
            f_7=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_7)*x));
            f_8=@(x)((exp(((-1)*gamma_g*(x.^2))-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_8)*x));
            
            Gauss(1,kk) = (1/pi)*integral(f_1,0,1000);
            Gauss(2,kk) = (1/pi)*integral(f_2,0,1000);
            Gauss(3,kk) = (1/pi)*integral(f_3,0,1000);
            Gauss(4,kk) = (1/pi)*integral(f_4,0,1000);
            Gauss_neg(1,kk) = (1/pi)*integral(f_5,0,1000);
            Gauss_neg(2,kk) = (1/pi)*integral(f_6,0,1000);
            Gauss_neg(3,kk) = (1/pi)*integral(f_7,0,1000);
            Gauss_neg(4,kk) = (1/pi)*integral(f_8,0,1000);
%             Astable(1,kk) = (1/pi)*integral(Sf_1,0,1000);
%             Astable(2,kk) = (1/pi)*integral(Sf_2,0,1000);
%             Astable(3,kk) = (1/pi)*integral(Sf_3,0,1000);
%             Astable(4,kk) = (1/pi)*integral(Sf_4,0,1000);
%             Astable_neg(1,kk) = (1/pi)*integral(Sf_5,0,1000);
%             Astable_neg(2,kk) = (1/pi)*integral(Sf_6,0,1000);
%             Astable_neg(3,kk) = (1/pi)*integral(Sf_7,0,1000);
%             Astable_neg(4,kk) = (1/pi)*integral(Sf_8,0,1000);
        end
        
        KJ=[-2.5:0.05:2.5];  
        samp_1=sum(Gauss(:,pos_gauss_start:pos_gauss_end));
        inv_samp_1=sum(Gauss(:,pos_alpha_right:end));
%         samp_2=sum(Astable);
        samp_3=sum(Gauss_neg(:,neg_gauss_start:neg_gauss_end));
        inv_samp_2=sum(Gauss_neg(:,neg_alpha_right:end));
%         samp_4=sum(Astable_neg);
        inv_samp_3=sum(Gauss(:,1:pos_alpha_left));
        inv_samp_4=sum(Gauss_neg(:,1:neg_alpha_left));
        model =  @(p,x) p(1)*(1/(sqrt(2*pi*p(3))))*exp(-(((x-p(2)).^2)/(2*(p(3)))));
        model_inv= @(g,y) g(1)./((g(2)*(abs(y).^g(3)))+g(4));
        cof_inv_1=nlinfit(KJ(pos_alpha_right:end), inv_samp_1, model_inv, [1 1 1 1]);
        cof_inv_2=nlinfit(KJ(neg_alpha_right:end), inv_samp_2, model_inv, [1 1 1 1]);
        cof_inv_3=nlinfit(KJ(1:pos_alpha_left), inv_samp_3, model_inv, [1 1 1 1]);
        cof_inv_4=nlinfit(KJ(1:neg_alpha_left), inv_samp_4, model_inv, [1 1 1 1]);
        coefEsts_1 = nlinfit(KJ(pos_gauss_start:pos_gauss_end), samp_1, model, [3 0.5 1]);
%         coefEsts_2 = nlinfit(KJ, samp_2, model_inv, [1 1 1 1]);
        coefEsts_3 = nlinfit(KJ(neg_gauss_start:neg_gauss_end), samp_3, model, [3 -0.5 1]);
%         coefEsts_4 = nlinfit(KJ, samp_4, model_inv, [3 1 1 1]);
        pos_inv=model_inv([cof_inv_1(1)    cof_inv_1(2)    cof_inv_1(3)  cof_inv_1(4)], KJ(pos_alpha_right:end));
        neg_inv=model_inv([cof_inv_2(1)    cof_inv_2(2)    cof_inv_2(3)  cof_inv_2(4)] , KJ(neg_alpha_right:end));
        pos_inv_2=model_inv([cof_inv_3(1)    cof_inv_3(2)    cof_inv_3(3)  cof_inv_3(4)], KJ(1:pos_alpha_left));
        neg_inv_2=model_inv([cof_inv_4(1)    cof_inv_4(2)    cof_inv_4(3)  cof_inv_4(4) ] , KJ(1:neg_alpha_left));
        Positive_1=model([coefEsts_1(1)    coefEsts_1(2)    coefEsts_1(3)], KJ(pos_gauss_start:pos_gauss_end));
%         Positive_2=model_inv([coefEsts_2(1)    coefEsts_2(2)    coefEsts_2(3) coefEsts_2(4)], KJ);
        negative_1=model([coefEsts_3(1)    coefEsts_3(2)    coefEsts_3(3)], KJ(neg_gauss_start:neg_gauss_end));
%         negative_2=model_inv([coefEsts_4(1)    coefEsts_4(2)    coefEsts_4(3)  coefEsts_4(4)], KJ);
        
        figure
        plot(KJ,sum(Gauss),'g+-');hold on;plot(KJ(pos_gauss_start:pos_gauss_end),Positive_1,'b+-');hold on;plot(KJ,sum(Gauss_neg),'g--');hold on ;plot(KJ(1:pos_alpha_left),(pos_inv_2),'k+-')
        hold on;plot(KJ(neg_gauss_start:neg_gauss_end),negative_1,'b--')
        hold on;plot(KJ(pos_alpha_right:end),(pos_inv),'r+-')
        hold on;plot(KJ(1:neg_alpha_left),(neg_inv_2),'k--');hold on;plot(KJ(neg_alpha_right:end),(neg_inv),'r--')
        legend('actual PDF for positive','lobe appr by gaussian ','actual PDF for negative','tail appr by alpha','entirely by apha','alpha estimated by gaussian')
        title(['SNR' '=' num2str(SNR)  '   '  '\kappa' '=' num2str(kappa)])
        

        LLR_Actual=log((sum(Gauss))./(sum(Gauss_neg)));
%         LLR_alpha=log(sum(Astable)./sum(Astable_neg));
        LLR_EST_B0=log([pos_inv_2  Positive_1  pos_inv] ./[neg_inv_2  negative_1  neg_inv]);       

        figure
        plot(KJ,LLR_Actual,'r*-');
        hold on;plot(KJ,LLR_EST_B0,'b-*');
%         hold on;plot(KJ,LLR_alpha,'gv--');
        legend('Exact LLR',' LLR by Gaussian and alpha Appr')
        xlabel('real(y)')
        ylabel('LLR(b_1)')
        title(['SNR' '=' num2str(SNR)  '   '  '\kappa' '=' num2str(kappa)])
        kappa_inx=kappa_inx+1;
    end
    SNR_inx=SNR_inx+1;
end

