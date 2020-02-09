function ZR = LLR_3(y,KJ,gamma_g,gamma_s,alpha,j)




full_pos = @(p,x) p(1)*(1/(sqrt(2*pi*p(3))))*(exp(-(((x-p(2)).^2)/(2*(p(3)))))+exp(-(((x+p(2)).^2)/(2*(p(3))))));
full_neg =  @(p,x) ((p(1)*(((1/(sqrt(2*pi*p(2))))*exp(-(((x).^2)/(2*(p(2))))))))+(p(3)*(((1/(sqrt(2*pi*p(4))))*exp(-(((x-p(5)).^2)/(2*(p(4))))))+((1/(sqrt(2*pi*p(4))))*exp(-(((x+p(5)).^2)/(2*(p(4)))))))));

% %  alpha=1.35 , snr=-2:2:6, kappa=1

% coff_pos_mat=[4.5608    0.6712    1.1613    0.8357    0.3027    2.0682    1.2028
%     4.7893    0.6006    0.7437    0.8132    0.1689    2.9819    1.2808
%     4.4720    0.6108    0.4621    0.9724    0.2524    3.4164    0.9182
%     5.2173    0.6090    0.2258    0.7114    0.1706    4.4525    1.1569
%     4.4506    0.6578    0.1770    1.1324    0.8079    3.6027    0.6762];
% coff_neg_mat=[8.9878    1.7375    0.7897    0.0461    3.5453    1.3458
%     9.6067    1.3666    0.7698    0.0548    4.3002    1.2843
%     6.1683    0.7731    1.2348    0.0253    5.5275    0.7738
%     1.8651    0.1426    1.5970    0.0229    7.0008    0.4811
%     1.7553    0.0895    1.6473    0.0147    8.2651    0.4812];

% SNR=-2:4:18  alpha=1.5   kappa=2
coff_pos_mat=[6.2236    0.6217    0.9814
    6.2444   -0.6182    0.3187
    6.2825    0.6182    0.1193
    6.3086    0.6189    0.0864
    6.4024    0.6200    0.0645
    6.4024    0.6200    0.0645
    6.5007    0.6209    0.0581];


coff_neg_mat=[ 0.0966    0.4860    6.2038    1.1882    0.6402
    5.2706    0.2739    3.6302    0.3224    1.0168
    6.3386    0.1212    3.0985    0.0902    1.0839
    6.3637    0.0882    3.0874    0.0556    1.0836
    6.4028    0.0645    3.1022    0.0311    1.0814
    6.4028    0.0645    3.1022    0.0311    1.0814
    6.4591    0.0570    3.1178    0.0311    1.0805];

coff_pos =coff_pos_mat(j,:);
coff_neg=coff_neg_mat(j,:);



F_1=full_pos([coff_pos(1) coff_pos(2) coff_pos(3)], KJ);
F_2=full_neg([coff_neg(1) coff_neg(2) coff_neg(3) coff_neg(4) coff_neg(5)], KJ);

Full_l=log(F_1./F_2);
% Threshold_left=-1.5;
% Threshold_rigth=1.5;
% half_KJ=((length(KJ))/2)+1;
%
% model_pos =  @(p,x) p(1)*(1/(sqrt(2*pi*p(3))))*(exp(-(((x-p(2)).^2)/(2*(p(3)))))+exp(-(((x+p(2)).^2)/(2*(p(3))))));
% model_neg =  @(p,x) (p(1)*(((1/(sqrt(2*pi*p(2))))*exp(-(((x).^2)/(2*(p(2))))))));
% model_inv= @(g,y) 1./((g(1)*(abs(y).^g(2)))+g(3));
% 
% 
% coefEsts_1 =[6.0482    0.6154    0.1589];
% coefEsts_3 =[ 12.5582    0.9540];
% 
% 
% 
% cof_inv_1 =[0.0936    6.4765    0.1689];
% cof_inv_2 =[ 0.0106    8.3811    0.2573];
% cof_inv_3=cof_inv_1;
% cof_inv_4=cof_inv_2;


% pos_inv=model_inv([cof_inv_1(1)    cof_inv_1(2)    cof_inv_1(3)  ], KJ(half_KJ+1:end));
% neg_inv=model_inv([cof_inv_2(1)    cof_inv_2(2)    cof_inv_2(3)  ] , KJ(half_KJ+1:end));
% pos_inv_2=model_inv([cof_inv_3(1)    cof_inv_3(2)    cof_inv_3(3) ], KJ(1:half_KJ));
% neg_inv_2=model_inv([cof_inv_4(1)    cof_inv_2(2)    cof_inv_4(3)  ] , KJ(1:half_KJ));
% Positive_1=model_pos([coefEsts_1(1)    coefEsts_1(2)    coefEsts_1(3)], KJ);
% negative_1=model_neg([coefEsts_3(1)    coefEsts_3(2)    coefEsts_3(3) ], KJ);
% 
% 
% LLR_EST_B0_right=log(pos_inv./neg_inv);
% LLR_EST_B0_center=log(Positive_1./negative_1);
% LLR_EST_B0_left=log(pos_inv_2./neg_inv_2);
% 
% KJ_LEFT=KJ(1:half_KJ);
% KJ_RIGHT=KJ(half_KJ+1:end);
% 
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
ZR=Full_l(c2);
        
end


