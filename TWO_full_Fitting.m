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
alpha=1.6;
SNR_inx=1;
KJ=[-2.5:0.1:2.5];
% %  KJ(1)=-2.5 ,  KJ(35)=-0.8  KJ(0)=51  KJ(67)=0.8  KJ(101)=2.5
% pos_gauss_start=41;
% pos_gauss_end=60;
% pos_alpha_left=40;
% pos_alpha_right=61;
%
% neg_gauss_start=3;
% neg_gauss_end=99;
% neg_alpha_left=2;
% neg_alpha_right=100;
% % 
% pos_gauss_start=22;
% pos_gauss_end=61;
% pos_alpha_left=pos_gauss_start-1;
% pos_alpha_right=pos_gauss_end+1;

neg_gauss_start=24;
neg_gauss_end=29;
neg_alpha_left=21;
neg_alpha_right=31;
KJ_half=(length(KJ)/2)+1;

pos_gauss_start=33;
pos_gauss_end=46;
pos_alpha_left=pos_gauss_start-1;
pos_alpha_right=pos_gauss_end+1;
% 
% neg_gauss_start=1;
% neg_gauss_end=21;
% neg_alpha_left=11;
% neg_alpha_right=21;

% pos_gauss_start=15;
% pos_gauss_end=36;
% pos_alpha_left=pos_gauss_start-1;
% pos_alpha_right=pos_gauss_end+1;
% 
% neg_gauss_start=15;
% neg_gauss_end=37;
% neg_alpha_left=15;
% neg_alpha_right=36;
% KJ_half=floor((length(KJ)/2))+1;
SNR_inx=0;
for SNR=15:5:15
    SNR_inx=SNR_inx+1;
    kappa_inx=1;
    for kappa=5:1:5
        gamma_s=((inv(10^(SNR/10)))*0.5)*(1/(kappa+1));
        gamma_g=kappa*gamma_s;
        gamma_s=0.1  ;
        gamma_g=0.1 ;
        
        Dist_1=0.1543 ;
        Dist_2=0.4629;
        Dist_3=-0.1543;
        Dist_4=-0.4629;
        
        
        Dist_5=0.7715;
        Dist_6=1.0801;
        Dist_7=-0.7715;
        Dist_8=-1.0801;
        
        kk=0;
        for KJ=[-2.5:0.1:2.5];
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
        
        KJ=[-2.5:0.1:2.5];
        samp_full=sum(Gauss(:,1:end));
        %         samp_full_neg=sum(Gauss_neg(:,1:half_KJ+intv));
        samp_full_neg=sum(Gauss_neg(:,1:end));
        samp_1=sum(Gauss(:,pos_gauss_start:pos_gauss_end));
        inv_samp_1=sum(Gauss(:,pos_alpha_right:end));
        samp_2=sum(Astable);
        samp_3=sum(Gauss_neg(:,neg_gauss_start:neg_gauss_end));
        inv_samp_2=sum(Gauss_neg(:,neg_alpha_right:end));
        samp_4=sum(Astable_neg);
        inv_samp_3=sum(Gauss(:,1:pos_alpha_left));
        inv_samp_4=sum(Gauss_neg(:,1:neg_alpha_left));
        model_pos =  @(p,x) p(1)*(1/(sqrt(2*pi*p(2))))*exp(-(((x).^2)/(2*(p(2)))));
        model_neg =  @(p,x) p(1)*(1/(sqrt(2*pi*p(3))))*(exp(-(((x-p(2)).^2)/(2*(p(3)))))+exp(-(((x+p(2)).^2)/(2*(p(3))))));
        %         model_neg =  @(p,x) 1*(((1/(sqrt(2*pi*p(1))))*(exp(-(((x-p(2)).^2)/(2*(p(1)))))+exp(-(((x+p(2)).^2)/(2*(p(1)))))))+p(3));
        model_inv= @(g,y) 1./((g(1)*(abs(y).^g(2)))+g(3));
        
        model_full_pos= @(h,z) ((h(1)*(1./((h(2)*(abs(z).^h(3)))+h(4)))) + (h(5)*(1/(sqrt(2*pi*h(6))))*exp(-(((z).^2)/(2*(h(6)))))));
        model_full_neg= @(h,w) ((h(1)*(1./((h(2)*(abs(w).^h(3)))+h(4)))) + (h(5)*(exp(-(((w-h(6)).^2)/(2*(h(7)))))+exp(-(((w+h(6)).^2)/(2*(h(7))))))));
        
        coff_full=nlinfit(KJ(:,1:end),samp_full, model_full_pos,[1 1 1 1 1 1 ]);
        vect_1(SNR_inx,:)=coff_full;
        coff_full_2=nlinfit(KJ(:,1:end),samp_full_neg, model_full_neg,[1 1 1 1 1 1 1]);
        %         coff_full_3=nlinfit(KJ,samp_full_neg_2, model_full_match,[1 1 1 1 1 1 1]);
        vect_2(SNR_inx,:)=coff_full_2;
        
        
        cof_inv_1=nlinfit(KJ(pos_alpha_right:end), inv_samp_1, model_inv, [1 1  1]);
        cof_inv_2=nlinfit(KJ(neg_alpha_right:end), inv_samp_2, model_inv, [1 1  1]);
        cof_inv_3=nlinfit(KJ(1:pos_alpha_left), inv_samp_3, model_inv, [1  1 1]);
        cof_inv_3=cof_inv_1;
        cof_inv_4=nlinfit(KJ(1:neg_alpha_left), inv_samp_4, model_inv, [1  1 1]);
        cof_inv_4=cof_inv_2;
        coefEsts_1 = nlinfit(KJ(pos_gauss_start:pos_gauss_end), samp_1, model_pos, [3  1]);
        coefEsts_2 = nlinfit(KJ, samp_2, model_inv, [1 1  1]);
        coefEsts_3 = nlinfit(KJ(neg_gauss_start:neg_gauss_end), samp_3, model_neg, [0.5 1 1]);
        coefEsts_4 = nlinfit(KJ, samp_4, model_inv, [3 1 1 ]);
        
        coff_full
        coff_full_2
        
        pos_full=model_full_pos([coff_full(1) coff_full(2) coff_full(3) coff_full(4) coff_full(5) coff_full(6) ], KJ);
        neg_full=model_full_neg([coff_full_2(1) coff_full_2(2) coff_full_2(3) coff_full_2(4) coff_full_2(5) coff_full_2(6) coff_full_2(7) ], KJ);
        
        pos_inv=model_inv([cof_inv_1(1)    cof_inv_1(2)    cof_inv_1(3)  ], KJ(KJ_half+1:end));
        neg_inv=model_inv([cof_inv_2(1)    cof_inv_2(2)    cof_inv_2(3)  ] , KJ(KJ_half+1:end));
        pos_inv_2=model_inv([cof_inv_3(1)    cof_inv_3(2)    cof_inv_3(3)  ], KJ(1:KJ_half));
        neg_inv_2=model_inv([cof_inv_4(1)    cof_inv_4(2)    cof_inv_4(3)   ] , KJ(1:KJ_half));
        Positive_1=model_pos([coefEsts_1(1)    coefEsts_1(2)    ], KJ);
        Positive_2=model_inv([coefEsts_2(1)    coefEsts_2(2)    coefEsts_2(3) ], KJ);
        negative_1=model_neg([coefEsts_3(1)    coefEsts_3(2)   coefEsts_3(3) ], KJ);
        negative_2=model_inv([coefEsts_4(1)    coefEsts_4(2)    coefEsts_4(3)  ], KJ);
        
        
        figure
        plot(KJ,sum(Gauss),'b-->');
        hold on;plot(KJ,Positive_1,'k--');
        hold on;plot(KJ, pos_full,'c-.+') ;
        hold on ; plot(KJ,sum(Gauss_neg),'r--<')
        hold on;plot(KJ,negative_1,'k-.')
        hold on ;plot(KJ, neg_full,'g-.+')
        hold on ;plot(zel(3:6), 0.35*ones(1,4),'b--o')
        hold on ;plot(zel(1:2), 0.35*ones(1,2),'r--o')
        hold on ;plot(zel(7:8), 0.35*ones(1,2),'r--o')
%         ylim([0.1  5.3])
        xt = [-1.1801  -0.8815  -0.5729  -0.2643  0.1643   0.4729   0.7815  1.1801];
        yt = 0.195*ones(1,4);
        yt_2 = 0.195*ones(1,2);
        str = {'-S_2','-S_1','S_1','S_2'}
        str_2={'S_4','S_3'};
        str_3={'-S_3','-S_4'};
        text(zel(7:8)-0.01,yt_2,str_2)
        hold on; text(zel(1:2)-0.02,yt_2,str_3)
        hold on; text(zel(3:6)-0.02,yt,str)
        xlabel('y_r')
        legend('Sum pdfs- $b_1$=1 (nom.)', 'Gaussian app. Eq.(10)', '$\hat{f}^{1}_1$', 'Sum pdfs-- $b_1$=0 (denom.)','Gaussian app. Eq.(12)', '$\hat{f}^{0}_1$','Location','northwest')
        set(gca, 'FontSize',14)
        set(gca, 'FontName', 'Times New Roman')
        set(legend,'Interpreter','latex')
        ylabel('sum pdfs')
        grid on
        
        LLR_Actual=log((sum(Gauss))./(sum(Gauss_neg)));
        LLR_EST_B0_left=log(pos_inv_2./neg_inv_2);
        LLR_EST_B0_center=log(Positive_1./negative_1);
        LLR_EST_B0_right=log(pos_inv./neg_inv);
        LLR_FULL=log(pos_full./neg_full);
        
        figure
        plot(KJ,LLR_Actual,'r-+');
        hold on;plot(KJ,(LLR_EST_B0_center),'k-.');
        hold on;plot(KJ, LLR_FULL, 'b-..') ;
        ylabel('\Lambda(b_1)')
        xlabel('y_r')
        legend('EXCAT LLR ', 'LLR App. by Gasuaain term',' LLR  App. by  $\hat{f}^{1}_1$  and   $\hat{f}^{0}_1 $','Location','northwest')
        set(gca, 'FontSize',13)
        set(gca, 'FontName', 'Times New Roman')
        set(legend,'Interpreter','latex')
%         ylim([-0.8  .7])
        grid on
        xlim([-2.6  2.6])
        %         figure
        %         plot(KJ,sum(Gauss),'gv--');hold on;plot(KJ,Positive_1,'b*-');hold on ;plot(KJ(1:KJ_half),(pos_inv_2),'k--')
        %         hold on ;plot(KJ, pos_full,'m-v')
        % %         hold on;plot(KJ(KJ_half+1:end),(pos_inv),'k--')
        %          legend('actual PDF f','lobe appr by gaussian ','tail','tail ')
        %         figure
        %         plot(KJ,sum(Gauss_neg),'rv--');hold on
        %         plot(KJ,negative_1,'b*-')
        %         hold on ;plot(KJ, neg_full,'m-v')
        %         hold on;plot(KJ(1:KJ_half),(neg_inv_2),'c-d');
        % %         hold on;plot(KJ(KJ_half+1:end),(neg_inv),'c-d')
        %         legend('actual PDF f','lobe appr by gaussian ','tail','tail ')
        %         grid on
        %         hold on;plot(KJ,samp_2,'k*-');hold on;plot(KJ,samp_4,'k*-')
        %         hold on;plot(KJ,Positive_2,'r*-');hold on;plot(KJ,negative_2,'r*-')
        
        %         title([num2str(SNR)  '--' num2str(kappa)])
        %         LLR_full=log(pos_full./neg_full);
        %
        %
        %         LLR_Actual=log((sum(Gauss))./(sum(Gauss_neg)));
        %         LLR_alpha=log(sum(Astable)./sum(Astable_neg));
        %         %         LLR_EST_B0=log([pos_inv_2  Positive_1  pos_inv] ./[neg_inv_2  negative_1  neg_inv]);
        %         LLR_EST_B0_left=log(pos_inv_2./neg_inv_2);
        %         LLR_EST_B0_center=log(Positive_1./negative_1);
        %         LLR_EST_B0_right=log(pos_inv./neg_inv);
        %
        %         figure
        %         plot(KJ,LLR_Actual,'r*-');
        %         hold on;plot(KJ,LLR_EST_B0_center,'b-*');
        %         hold on;plot( KJ(1:KJ_half),LLR_EST_B0_left,'g-d');
        %         hold on;plot(KJ(KJ_half+1:end),LLR_EST_B0_right,'g-v');
        %         hold on; plot(KJ,LLR_full,'c-.s')
        %         legend('Exact LLR','Gaussian',' tail left',' tail right')
        %         xlabel('real(y)')
        %         ylabel('LLR(b_1)')
        %         grid on
%         title([num2str(SNR)  '--' num2str(kappa)])
        kappa_inx=kappa_inx+1;
    end
end

