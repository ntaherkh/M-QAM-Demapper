clc
clear
close all
KK=[-1.0801    -0.7715  -0.4629  -0.1543  0.1543 0.4629  0.7715  1.0801];

delta=0;
beta_skew=0;
zel=[-1.0801  -0.7715  -0.4629 -0.1543  0.1543  0.4629  0.7715  1.0801];
sefr=zeros(1,8);

SNR_inx=1;
KJ=[-2.5:0.1:2.5];

pos_alpha_left=21;% 37;
pos_gauss_end=27;% 64;
pos_gauss_start=pos_alpha_left+1;% 38;
pos_alpha_right=pos_gauss_end+1;% 65;
neg_alpha_right=35;
neg_alpha_left=23;
neg_gauss_start=neg_alpha_left+1;
neg_gauss_end=neg_alpha_right-1;


% pos_alpha_left=11;% 37;
% pos_gauss_end=33;% 64;
% pos_gauss_start=pos_alpha_left+1;% 38;
% pos_alpha_right=pos_gauss_end+1;% 65;

% pos_alpha_left=11;% 37;
% pos_gauss_end=33;% 64;
% pos_gauss_start=pos_alpha_left+1;% 38;
% pos_alpha_right=pos_gauss_end+1;% 65;
% 


half_KJ=floor(0.5*(length(KJ)))+1;
intv=25;

% alpha_vect=[1.1:0.2:1.9];
alpha_inx=0;
for alpha=1.6:0.1:1.6
    alpha_inx=alpha_inx+1;
    SNR_inx=0;
    for SNR=15:1:15
        SNR_inx=SNR_inx+1;
        kappa_inx=0;
        for kappa=5:1:5
            kappa_inx=kappa_inx+1;
%             gamma_s=((inv(10^(SNR/10)))*0.5)*(1/(kappa+1));
%             gamma_g=kappa*gamma_s;
            gamma_s=0.1;
            gamma_g=0.1;
            Dist_1=0.1543 ;
            Dist_2=0.4629;
            Dist_3=0.7715;
            Dist_4=1.0801;
            
            Dist_5=-0.1543 ;
            Dist_6=-0.4629 ;
            Dist_7=-0.7715;
            Dist_8=-1.0801;
            kk=0;
            for KJ=[-2.5:0.1:2.5];
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
%             gamma_g=0.148;
%             gamma_s=0.11;
            
            KJ=[-2.5:0.1:2.5];
            samp_full=sum(Gauss(:,half_KJ-intv:end));
            samp_full_neg=sum(Gauss_neg(:,1:half_KJ+intv));
            samp_full_neg_2=sum(Gauss_neg(:,1:end));
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
            model_full_match= @(h,z) ((h(1)*(1./((h(2)*(abs(z).^h(3)))+h(4)))) + (h(5)*(1/(sqrt(2*pi*h(7))))*exp(-(((z-h(6)).^2)/(2*(h(7)))))));
            model_full_2= @(h,z) ((h(1)*(1/(sqrt(2*pi*h(3))))*exp(-(((z+h(2)).^2)/(2*(h(3)))))));
            
            coff_full=nlinfit(KJ(:,half_KJ-intv:end),samp_full, model_full_match,[1 1 1 1 1 1 1]);
%             coff_full_3=nlinfit(KJ(:,half_KJ-intv:end),samp_full, model_full_2,[1 1 1 1 1 1 1]);
%             vect_1(SNR_inx,:)=
            kappa 
            alpha
            coff_full
%             coff_full_2=nlinfit(KJ(:,1:half_KJ+intv),samp_full_neg, model_full_match,[1 1 1 1 1 1 1]);
            coff_full_3=nlinfit(KJ,samp_full, model_full_2,[1 1 1 ]);
            
            coff_full_2=coff_full;
            coff_full_2(6)=-coff_full(6);
            vect_2(SNR_inx,:)=coff_full_2;
            cof_inv_1=nlinfit(KJ(pos_alpha_right:end), inv_samp_1, model_inv, [1 1 1 1]);
            cof_inv_2=nlinfit(KJ(neg_alpha_right:end), inv_samp_2, model_inv, [1 1 1 1]);
            cof_inv_3=nlinfit(KJ(1:pos_alpha_left), inv_samp_3, model_inv, [1 1 1 1]);
            cof_inv_4=nlinfit(KJ(1:neg_alpha_left), inv_samp_4, model_inv, [1 1 1 1]);
            cof_inv_4=cof_inv_1;
            cof_inv_3=cof_inv_2;
            coefEsts_1 = nlinfit(KJ(pos_gauss_start:pos_gauss_end), samp_1, model, [3 0.5 1]);
%             coefEsts_2 = nlinfit(KJ, samp_2, model_inv, [1 1 1 1]);
            coefEsts_3 = nlinfit(KJ(neg_gauss_start:neg_gauss_end), samp_3, model, [3 -0.5 1]);
%                     coefEsts_4 = nlinfit(KJ, samp_4, model_inv, [3 1 1 1]);
            pos_full=model_full_match([coff_full(1) coff_full(2) coff_full(3) coff_full(4) coff_full(5) coff_full(6) coff_full(7)], KJ);
            pos_full_2=model_full_2([coff_full_3(1) coff_full_3(2) coff_full_3(3)], KJ);
            neg_full=model_full_match([coff_full_2(1) coff_full_2(2) coff_full_2(3) coff_full_2(4) coff_full_2(5) coff_full_2(6) coff_full_2(7)], KJ);
%             neg_full_2=model_full_match([coff_full_3(1) coff_full_3(2) coff_full_3(3) coff_full_3(4) coff_full_3(5) coff_full_3(6) coff_full_3(7)], KJ);
            pos_inv=model_inv([cof_inv_1(1)    cof_inv_1(2)    cof_inv_1(3)  cof_inv_1(4) ], KJ(half_KJ+1:end));
            neg_inv=model_inv([cof_inv_2(1)    cof_inv_2(2)    cof_inv_2(3)  cof_inv_2(4)] , KJ(half_KJ+1:end));
            pos_inv_2=model_inv([cof_inv_3(1)    cof_inv_3(2)    cof_inv_3(3)  cof_inv_3(4)], KJ(1:half_KJ));
            neg_inv_2=model_inv([cof_inv_4(1)    cof_inv_4(2)    cof_inv_4(3)  cof_inv_4(4) ] , KJ(1:half_KJ));
            Positive_1=model([coefEsts_1(1)    coefEsts_1(2)    coefEsts_1(3)], KJ);
            %         Positive_2=model_inv([coefEsts_2(1)    coefEsts_2(2)    coefEsts_2(3) coefEsts_2(4)], KJ);
            negative_1=model([coefEsts_3(1)    coefEsts_3(2)    coefEsts_3(3)], KJ);
            %         negative_2=model_inv([coefEsts_4(1)    coefEsts_4(2)    coefEsts_4(3)  coefEsts_4(4)], KJ);
            
            tensor_cof_inv_1(alpha_inx,kappa_inx,:)=cof_inv_1;
            tensor_cof_inv_2(alpha_inx,kappa_inx,:)=cof_inv_2;
            tensor_cof_inv_3(alpha_inx,kappa_inx,:)=cof_inv_3;
            tensor_cof_inv_4(alpha_inx,kappa_inx,:)=cof_inv_4;
            tensor_coefEsts_1(alpha_inx,kappa_inx,:)=coefEsts_1;
            tensor_coefEsts_2(alpha_inx,kappa_inx,:)=coefEsts_3;
            
            figure
            plot(KJ,sum(Gauss),'b-->');
            hold on;plot(KJ,Positive_1,'k--');
            hold on;plot(KJ, pos_full,'c-.+') ;
%             hold on;plot(KJ, pos_full_2,'c--');
            hold on ; plot(KJ,sum(Gauss_neg),'r--<')
            hold on;plot(KJ,negative_1,'k-.')
            hold on ;plot(KJ, neg_full,'g-.+')
            hold on ;plot(zel(1:4), 0.22*ones(1,4),'r--o')
            hold on ;plot(zel(5:8), 0.22*ones(1,4),'b--o')
            xt = [-1.1801  -0.8815  -0.5729  -0.2643  0.1643   0.4729   0.7815  1.1801];
            yt = 0.165*ones(1,4);
            yt_2 = 0.165*ones(1,4);
            str = {'-S_4','-S_3','-S_2','-S_1'};
            str_2={'S_1','S_2','S_3','S_4'};
            text(zel(1:4)-0.05,yt,str)
            hold on; text(zel(5:8)-0.02,yt_2,str_2)
            xlabel('y_r')
            legend('Sum pdfs- $b_0$=1 (nom.)', 'Gaussian app. Eq.(8)', '$\hat{f}^{1}_0$', 'Sum pdfs-- $b_0$=0 (denom.)','Gaussian app. Eq.(8)', '$\hat{f}^{0}_0$','Location','northwest')
            set(gca, 'FontSize',13)
            set(gca, 'FontName', 'Times New Roman')
            set(legend,'Interpreter','latex')
            ylabel('sum pdfs')
            xlim([-2.5  2.5])
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
            ylabel('\Lambda(b_0)')
            xlabel('y_r')
            legend('EXCAT LLR ', 'LLR App. by Gasuaain term',' LLR  App. by  $\hat{f}^{1}_0$  and   $\hat{f}^{0}_0 $','Location','northwest')
            set(gca, 'FontSize',14)
            set(gca, 'FontName', 'Times New Roman')
            set(legend,'Interpreter','latex')
            xlim([-2.6  2.6])
            
            
            grid on
            
%             hold on ;plot(zel, 0.2*ones(1,length(zel)),'k-o')
            xt = [-1.1801  -0.8815  -0.5729  -0.2643  0.1643   0.4729   0.7815  1.1801];
%             yt = 0.175*ones(1,8);
%             str = {'-S_4','-S_3','-S_2','-S_1','S_1','S_2','S_3','S_4'};
% %             text(zel,yt,str)
%             xlabel('y_r')
%             legend('Sum pdfs- $b_0$=1 (nom.)', 'Gaussian app. Eq.(8)', '$\hat{f}^{1}_0$', 'Sum pdfs-- $b_0$=0 (denom.)','Gaussian app. Eq.(8)', '$\hat{f}^{0}_0$')
%             set(gca, 'FontSize',14)
%             set(gca, 'FontName', 'Times New Roman')
%             set(legend,'Interpreter','latex')
%             ylim([0.1  1.7])
            
            
            %
%             hold on;plot(KJ(1:pos_alpha_left),(pos_inv_2(1:pos_alpha_left)),'r--')
%             hold on;plot(KJ(pos_alpha_right:end),(pos_inv(end+1-length(KJ(pos_alpha_right:end)):end)),'r--')
%             legend('f(u)','$\hat{f}^{1}_0$')
%            
            grid on
%             figure
%             plot(KJ,sum(Gauss_neg),'g-*');hold on
% 
%             hold on;plot(KJ,negative_1,'b--')
%             hold on;plot(KJ(1:neg_alpha_left),(neg_inv_2(1:length(KJ(1:neg_alpha_left)))),'r--');hold on;plot(KJ(neg_alpha_right:end),(neg_inv(end-(length(KJ(neg_alpha_right:end)))+1:end)),'r--')
%             legend('actual PDF for negative','center appr by gaussian ','left tail appr by alpha','right tail by alpha')
%             title(['\alpha' '=' num2str(alpha)  '   ' 'SNR' '=' num2str(SNR)  '   '  '\kappa' '=' num2str(kappa)])
%             grid on
%             LLR_Actual=log((sum(Gauss))./(sum(Gauss_neg)));
%             %         LLR_alpha=log(sum(Astable)./sum(Astable_neg));
% %             LLR_EST_B0=log([pos_inv_2  Positive_1  pos_inv] ./[neg_inv_2  negative_1  neg_inv]);
%             LLR_EST_B0_left=log(pos_inv_2./neg_inv_2);
%             LLR_EST_B0_center=log(Positive_1./negative_1);
%             LLR_EST_B0_right=log(pos_inv./neg_inv);
%             LLR_FULL=log(pos_full./neg_full);
%             figure
%             plot(KJ,LLR_Actual,'r*-');
%             hold on;plot(KJ,(LLR_EST_B0_center),'g+-')
%             hold on;plot(KJ(1:half_KJ),(LLR_EST_B0_left),'k-d')
%             hold on;plot(KJ(half_KJ+1:end),(LLR_EST_B0_right),'k-v')
%             hold on; plot(KJ, LLR_FULL,'m-.*')
%             legend('Exact LLR',' LLR by Gaussian , LEFT tail' ,'RIGHT tail')
%             xlabel('real(y)')
%             ylabel('LLR(b_1)')
%             title(['\alpha' '=' num2str(alpha)  '   ' 'SNR' '=' num2str(SNR)  '   '  '\kappa' '=' num2str(kappa)])
%             grid on
% %                         plot(KJ,LLR_EST_B0,'b-*');
% %                     hold on;plot(KJ,LLR_alpha,'gv--');
% %                         for i=1:half_KJ
% %                             LLR_ESTMEAN_B0(i)=max(LLR_EST_B0_left(i),LLR_EST_B0_center(i));
% %                         end
% %                         for i=half_KJ+1:1:length(KJ)
% %                             LLR_ESTMEAN_B0(i)=min(LLR_EST_B0_right(i-half_KJ),LLR_EST_B0_center(i));
% %                         end
% %                         expected_error(alpha_inx,kappa_inx)=sum((LLR_Actual-LLR_ESTMEAN_B0).^2)/length(LLR_Actual);
%             [cof_inv_1    cof_inv_2] 
%             [cof_inv_3    cof_inv_4]
%             [coefEsts_1   coefEsts_3]
        end
    end
end
% expected_error

