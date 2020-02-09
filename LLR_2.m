function ZR = LLR_2(y,KJ,gamma_g,gamma_s,alpha,j)



Threshold_left=-1.5;
Threshold_rigth=1.5;

model_full_pos= @(h,z) ((h(1)*(1./((h(2)*(abs(z).^h(3)))+h(4)))) + (h(5)*(1/(sqrt(2*pi*h(6))))*exp(-(((z).^2)/(2*(h(6)))))));
model_full_neg= @(h,w) ((h(1)*(1./((h(2)*(abs(w).^h(3)))+h(4)))) + (h(5)*(exp(-(((w-h(6)).^2)/(2*(h(7)))))+exp(-(((w+h(6)).^2)/(2*(h(7))))))));

%  alpha=1.35 , snr=-2:2:6, kappa=1
% coff_full_2_mat=[0.8356    0.2433    2.0684    1.3477    1.6988    0.9547    1.1487
%      0.9625    0.1098    3.3824    1.0230    1.9442    0.9397    0.7295
%      1.2182    0.2047    3.4796    0.8321    2.3215    1.0073    0.4188
%      0.8457    0.0341    5.4471    1.2847    3.7402    0.9214    0.2707
%       0.5956    0.0598    5.1286    1.7078    5.1257    0.9250    0.1859];
%
%
% coff_full_mat=[0.8109    0.3436    2.0617    1.1587    9.2703    1.3540
%     0.5178    0.0254    4.2502    1.9036   10.7793    0.8537
%     0.7379    0.2687    3.2593    1.0675   10.0384    0.5444
%     0.9808    0.5016    3.7811    0.7062    9.0526    0.3648
%     1.2037    0.9129    4.4070    0.5052    7.9284    0.2653];

%  KAPPA=5 SNR=[-2:2:8] ALPHA=1.7
% coff_full_mat=[        0.5691    0.5726    2.0697    2.0902    1.8663    0.9537    1.4723
%     0.6494    0.0156    4.9911    1.6630    2.2037    0.8894    0.9042
%     0.7003    0.0215    5.5590    1.7131    2.7904    0.9070    0.5849
%     1.3166    0.3774    3.9818    0.7158    2.7660    1.0846    0.3196
%     0.1850    0.4241    3.1302    5.3087    4.8414    0.9254    0.2611
%     0.1850    0.4241    3.1302    5.3087    4.8414    0.9254    0.2611];
% 
% 
% coff_full_2_mat=[    0.6288    0.4624    2.0578    1.8650    1.8317    0.9557    1.4163
%     0.6697    0.0176    4.8344    1.6123    2.2067    0.8913    0.8679
%     0.7279    0.0245    5.3639    1.6377    2.8012    0.9080    0.5585
%     1.1870    0.4224    3.8761    0.8270    3.0500    1.0344    0.3295
%     0.2007    0.2696    3.2499    3.2311    4.9028    0.9254    0.2491
%     0.0068    386421    -31.2672  0.4484    6.0840    0.9253    0.1684];

%     snr=-2:4:18 alpha=1.5 kappa=2    


coff_full_mat=[    0.3309    0.2022    3.3546    3.0053   12.0022    1.0581
     0.6058    0.9419    5.0853    1.1439   11.3708    0.4141
    2.1380    0.4625    3.1219    0.1578   -9.7847    0.9676
    2.1085    0.7618    3.8713    0.1532   -7.6397    0.6128
    1.7393    4.7522    7.1539    0.1912    0.7770    0.1441
     1.7393    4.7522    7.1539    0.1912    0.7770    0.1441
    1.7136   13.9582    9.2345    0.1757    0.2413    0.1121];
    


coff_full_2_mat=[0.6548    0.0162    4.9377    1.6477    2.2060    0.8901    0.8928
    0.3323    0.1140    4.8204    3.3010    4.3114    0.9249    0.3147
   -0.1620   -2.0000   -1.0000    1.9465    7.3000    0.9256    0.1224
    -0.5810    0.1000   -3.0000    5.1470    8.6444    0.9217    0.0875
     -0.2852  -10.0000   -5.0000   16.3340   10.0673    0.9262    0.0642
     -0.2852  -10.0000   -5.0000   16.3340   10.0673    0.9262    0.0642
     -0.2852   -7.0000   -5.0000   16.3340   10.0673    0.9262    0.0642];

coff_full_2 =coff_full_2_mat(j,:);
coff_full =coff_full_mat(j,:);

pos_full=model_full_pos([coff_full(1) coff_full(2) coff_full(3) coff_full(4) coff_full(5) coff_full(6) ], KJ);
neg_full=model_full_neg([coff_full_2(1) coff_full_2(2) coff_full_2(3) coff_full_2(4) coff_full_2(5) coff_full_2(6) coff_full_2(7) ], KJ);


LLR_FULL=real(log(pos_full./neg_full));
% KJ_half=((length(KJ))/2)+1;
%
% model_pos =  @(p,x) p(1)*(1/(sqrt(2*pi*p(2))))*exp(-(((x).^2)/(2*(p(2)))));
% model_neg =  @(p,x) p(1)*(1/(sqrt(2*pi*p(3))))*(exp(-(((x-p(2)).^2)/(2*(p(3)))))+exp(-(((x+p(2)).^2)/(2*(p(3))))));
% model_inv= @(g,y) g(1)./((g(2)*(abs(y).^g(3)))+g(4));
%
% coefEsts_1 =[12.3864    0.3029];
% coefEsts_3 =[6.1111    0.9215    0.1692];
%
%
% cof_inv_1 =[0.5201    4.7378    0.1489];
% cof_inv_3=cof_inv_1;
%
% cof_inv_2 = [0.0046    9.9563    0.1859];
% cof_inv_4=cof_inv_2;
%
%
% pos_inv=model_inv([cof_inv_1(1)    cof_inv_1(2)    cof_inv_1(3)  cof_inv_1(4)], KJ(KJ_half+1:end));
% neg_inv=model_inv([cof_inv_2(1)    cof_inv_2(2)    cof_inv_2(3)  cof_inv_2(4)] , KJ(KJ_half+1:end));
% pos_inv_2=model_inv([cof_inv_3(1)    cof_inv_3(2)    cof_inv_3(3)  cof_inv_3(4)], KJ(1:KJ_half));
% neg_inv_2=model_inv([cof_inv_4(1)    cof_inv_4(2)    cof_inv_4(3)  cof_inv_4(4) ] , KJ(1:KJ_half));
% Positive_1=model_pos([coefEsts_1(1)    coefEsts_1(2)    ], KJ);
% negative_1=model_neg([coefEsts_3(1)    coefEsts_3(2)   coefEsts_3(3) ], KJ);
%
% LLR_EST_B0_right=log(pos_inv./neg_inv);
% LLR_EST_B0_center=log(Positive_1./negative_1);
% LLR_EST_B0_left=log(pos_inv_2./neg_inv_2);
%
% KJ_LEFT=KJ(1:KJ_half);
% KJ_RIGHT=KJ(KJ_half+1:end);

% if ((y >= Threshold_left)&&(y <=Threshold_rigth))
%     [c1,c2]=min(abs(KJ-y));
%     ZR=LLR_EST_B0_center(c2)
% elseif( y<Threshold_left)
%     [d1,d2]=min(abs(KJ_LEFT-y));
%     ZR=LLR_EST_B0_left(d2)
% elseif ( y>Threshold_rigth)
%     [e1,e2]=min(abs(KJ_RIGHT-y));
%     ZR=LLR_EST_B0_right(e2)
%
% end
[~,c2]=min(abs(KJ-y));
ZR=LLR_FULL(c2);

end


