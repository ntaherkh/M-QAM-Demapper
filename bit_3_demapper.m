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
pos_gauss_start=32;
pos_gauss_end=70;
pos_alpha_left=31;
pos_alpha_right=71;

neg_gauss_start=27;
neg_gauss_end=75;
neg_alpha_left=26;
neg_alpha_right=76;

for SNR=5:5:15
    kappa_inx=1;
    for kappa=1:2:5
        gamma_g=(inv(10^(SNR/10)))*(1/(kappa+1));
        gamma_s=kappa*gamma_g;
        
        Dist_1=-0.7715;
        Dist_2=-0.4629;
        Dist_3=0.46291;
        Dist_4=0.7715;
       
        Dist_5=-1.0801;
        Dist_6=-0.1543;
        Dist_7=0.1543;
        Dist_8=1.0801;
        
        kk=0;
        for KJ=[-2.5:0.05:2.5];
            kk=kk+1;
            Sf_1=@(x)((exp(-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_1)*x));
            Sf_2=@(x)((exp(-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_2)*x));
            Sf_3=@(x)((exp(-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_3)*x));
            Sf_4=@(x)((exp(-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_4)*x));
            Sf_5=@(x)((exp(-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_5)*x));
            Sf_6=@(x)((exp(-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_6)*x));
            Sf_7=@(x)((exp(-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_7)*x));
            Sf_8=@(x)((exp(-gamma_s*(abs(x).^alpha))).*cos((KJ-Dist_8)*x));
            
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
            Astable(1,kk) = integral(Sf_1,0,1000);
            Astable(2,kk) = integral(Sf_2,0,1000);
            Astable(3,kk) = integral(Sf_3,0,1000);
            Astable(4,kk) = integral(Sf_4,0,1000);
            Astable_neg(1,kk) = integral(Sf_5,0,1000);
            Astable_neg(2,kk) = integral(Sf_6,0,1000);
            Astable_neg(3,kk) = integral(Sf_7,0,1000);
            Astable_neg(4,kk) = integral(Sf_8,0,1000);
        end
        
        KJ=[-2.5:0.05:2.5];  
        samp_1=sum(Gauss(:,pos_gauss_start:pos_gauss_end));
        inv_samp_1=sum(Gauss(:,pos_alpha_right:end));
        samp_2=sum(Astable);
        samp_3=sum(Gauss_neg(:,neg_gauss_start:neg_gauss_end));
        inv_samp_2=sum(Gauss_neg(:,neg_alpha_right:end));
        samp_4=sum(Astable_neg);
        inv_samp_3=sum(Gauss(:,1:pos_alpha_left));
        inv_samp_4=sum(Gauss_neg(:,1:neg_alpha_left));
        model_pos =  @(p,x) p(1)*(1/(sqrt(2*pi*p(3))))*(exp(-(((x-p(2)).^2)/(2*(p(3)))))+exp(-(((x+p(2)).^2)/(2*(p(3))))));
        model_neg =  @(p,x) p(1)*(((1/(sqrt(2*pi*p(2))))*exp(-(((x).^2)/(2*(p(2))))))+((1/(sqrt(2*pi*p(3))))*(exp(-(((x-p(4)).^2)/(2*(p(3)))))+exp(-(((x+p(4)).^2)/(2*(p(3))))))));
        model_inv= @(g,y) g(1)./((g(2)*(abs(y).^g(3)))+g(4));
        cof_inv_1=nlinfit(KJ(pos_alpha_right:end), inv_samp_1, model_inv, [1 1 1 1]);
        cof_inv_2=nlinfit(KJ(neg_alpha_right:end), inv_samp_2, model_inv, [1 1 1 1]);
        cof_inv_3=nlinfit(KJ(1:pos_alpha_left), inv_samp_3, model_inv, [1 1 1 1]);
        cof_inv_4=nlinfit(KJ(1:neg_alpha_left), inv_samp_4, model_inv, [1 1 1 1]);
        coefEsts_1 = nlinfit(KJ(pos_gauss_start:pos_gauss_end), samp_1, model_pos, [3 0.5 1]);
        coefEsts_2 = nlinfit(KJ, samp_2, model_inv, [1 1 1 1]);
        coefEsts_3 = nlinfit(KJ(neg_gauss_start:neg_gauss_end), samp_3, model_neg, [3 0.5 0.5 1]);
        coefEsts_4 = nlinfit(KJ, samp_4, model_inv, [3 1 1 1]);
        pos_inv=model_inv([cof_inv_1(1)    cof_inv_1(2)    cof_inv_1(3)  cof_inv_1(4)], KJ(pos_alpha_right:end));
        neg_inv=model_inv([cof_inv_2(1)    cof_inv_2(2)    cof_inv_2(3)  cof_inv_2(4)] , KJ(neg_alpha_right:end));
        pos_inv_2=model_inv([cof_inv_3(1)    cof_inv_3(2)    cof_inv_3(3)  cof_inv_3(4)], KJ(1:pos_alpha_left));
        neg_inv_2=model_inv([cof_inv_4(1)    cof_inv_4(2)    cof_inv_4(3)  cof_inv_4(4) ] , KJ(1:neg_alpha_left));
        Positive_1=model_pos([coefEsts_1(1)    coefEsts_1(2)    coefEsts_1(3)], KJ(pos_gauss_start:pos_gauss_end));
        Positive_2=model_inv([coefEsts_2(1)    coefEsts_2(2)    coefEsts_2(3) coefEsts_2(4)], KJ);
        negative_1=model_neg([coefEsts_3(1)    coefEsts_3(2)    coefEsts_3(3) coefEsts_3(4)], KJ(neg_gauss_start:neg_gauss_end));
        negative_2=model_inv([coefEsts_4(1)    coefEsts_4(2)    coefEsts_4(3)  coefEsts_4(4)], KJ);
        
        figure
        plot(KJ,sum(Gauss),'gv--');hold on;plot(KJ(pos_gauss_start:pos_gauss_end),Positive_1,'b*-');hold on;plot(KJ,sum(Gauss_neg),'rv--');hold on ;plot(KJ(1:pos_alpha_left),(pos_inv_2),'k--')
        hold on;plot(KJ(neg_gauss_start:neg_gauss_end),negative_1,'b*-')
        hold on;plot(KJ(pos_alpha_right:end),(pos_inv),'k--')
        hold on;plot(KJ(1:neg_alpha_left),(neg_inv_2),'k--');hold on;plot(KJ(neg_alpha_right:end),(neg_inv),'k--')
%         hold on;plot(KJ,samp_2,'k*-');hold on;plot(KJ,samp_4,'k*-')
%         hold on;plot(KJ,Positive_2,'r*-');hold on;plot(KJ,negative_2,'r*-')
        legend('actual PDF for positive','lobe appr by gaussian ','actual PDF for negative','tail appr by alpha','entirely by apha','alpha estimated by gaussian')
        title([num2str(SNR)  '--' num2str(kappa)])
        

        LLR_Actual=log((sum(Gauss))./(sum(Gauss_neg)));
        LLR_alpha=log(sum(Astable)./sum(Astable_neg));
        LLR_EST_B0=log([pos_inv_2  Positive_1  pos_inv] ./[neg_inv_2  negative_1  neg_inv]);       

        figure
        plot(KJ,LLR_Actual,'r*-');
        hold on;plot(KJ,LLR_EST_B0,'b-*');
        hold on;plot(KJ,LLR_alpha,'gv--');
        legend('Exact LLR','Exact only Alpha',' LLR by Gaussian Appr','LLR by tail Appr')
        xlabel('real(y)')
        ylabel('LLR(b_1)')
        title([num2str(SNR)  '--' num2str(kappa)])
        kappa_inx=kappa_inx+1;
    end
    SNR_inx=SNR_inx+1;
end

