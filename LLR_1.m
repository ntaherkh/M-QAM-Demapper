function ZR = LLR_1(y,KJ,gamma_g,gamma_s,alpha,j)



Threshold_left=-1;
Threshold_rigth=1;
%  alpha=1.35 , snr=-2:2:6, kappa=1
% coff_full_mat =[0.3556    0.0286    3.4167    6.4418    3.4593    0.6187    1.3884
%     0.2882    0.1517    2.4761    5.8476    3.5932    0.6196    0.8944
%     0.2524    0.2594    2.4255    6.5202    3.7122    0.6188    0.6029
%     0.1557    2.6768    1.2062    1.6175    3.7394    0.6300    0.4207
%     0.2604    9.7984    2.4393    0.9745    3.6771    0.6561    0.3013];

% SNR=-2:2:8  alpha=1.7    kappa=5


% coff_full_mat=[
%     0.2194    1.9532    2.1219    7.9433    3.8605    0.6269    1.6533
%     0.1823    3.2220    2.1395    7.0213    3.8966    0.6249    1.0859
%     0.1903    7.3647    2.2853    6.5085    3.9201    0.6250    0.7308
%     2.0752   75.8207    2.8318   21.5601    3.8470    0.6397    0.5006
%     1.0848   42.5450    3.1530    4.2123    3.7106    0.6647    0.3474
%     1.7226    6426      4.3100   73.7479   3.7833     0.6580    0.2683];


% SNR= [0 5 10 12 14.5 15 16]
%     snr=-2:4:18 alpha=1.5 kappa=2
coff_full_mat=[
    0.2191    2.8241    2.1078    6.9005    3.8653    0.6259    1.0725
    1.1756   35.3591    2.9128    6.3467    3.7503    0.6557    0.4103
    0.2519    4.0000    6.0237    0.7335    3.8326    0.6538    0.2219
    0.1651   25.0000    6.6126    0.4444    3.8926    0.6500    0.1996
   0.1341   60.0000    6.7402    0.3584    3.9794    0.6446    0.1884
   0.1341   60.0000    6.7402    0.3584    3.9794    0.6446    0.1884
   0.1341   60.0000    6.7402    0.3584    3.9794    0.6446    0.1884];



coff_full=coff_full_mat(j,:);
coff_full_2=coff_full;
coff_full_2(6)=-coff_full(6);
model_full_match= @(h,z) ((h(1)*(1./((h(2)*(abs(z).^h(3)))+h(4)))) + (h(5)*(1/(sqrt(2*pi*h(7))))*exp(-(((z-h(6)).^2)/(2*(h(7)))))));
pos_full=model_full_match([coff_full(1) coff_full(2) coff_full(3) coff_full(4) coff_full(5) coff_full(6) coff_full(7)], KJ);
neg_full=model_full_match([coff_full_2(1) coff_full_2(2) coff_full_2(3) coff_full_2(4) coff_full_2(5) coff_full_2(6) coff_full_2(7)], KJ);
LLR_FULL=log(pos_full./neg_full);

% half_KJ=((length(KJ))/2)+1;
%
% model_inv= @(g,y) g(1)./((g(2)*(abs(y).^g(3)))+g(4));
% model =  @(p,x) p(1)*(1/(sqrt(2*pi*p(3))))*exp(-(((x-p(2)).^2)/(2*(p(3)))));
%
% coefEsts_1 =[3.9666    0.6226    0.3071];
% coefEsts_3=coefEsts_1;
% coefEsts_3(2)=(-1)*coefEsts_1(2);
%
% cof_inv_1 =[ 0.0403    7.7044    0.4388];
% cof_inv_2 =[ 21.9261    1.9138   -4.4843];
% cof_inv_3=cof_inv_2;
% cof_inv_4=cof_inv_1;
%
% coff_full_2=nlinfit(KJ(:,1:half_KJ+intv),samp_full_neg, model_full_match,[1 1 1 1 1 1 1]);



% Positive_1=model([coefEsts_1(1)    coefEsts_1(2)    coefEsts_1(3)], KJ);
% negative_1=model([coefEsts_3(1)    coefEsts_3(2)    coefEsts_3(3)], KJ);
% pos_inv=model_inv([cof_inv_1(1)    cof_inv_1(2)    cof_inv_1(3)  cof_inv_1(4) ], KJ(half_KJ+1:end));
% neg_inv=model_inv([cof_inv_2(1)    cof_inv_2(2)    cof_inv_2(3)  cof_inv_2(4)] , KJ(half_KJ+1:end));
% pos_inv_2=model_inv([cof_inv_3(1)    cof_inv_3(2)    cof_inv_3(3)  cof_inv_3(4)], KJ(1:half_KJ));
% neg_inv_2=model_inv([cof_inv_4(1)    cof_inv_4(2)    cof_inv_4(3)  cof_inv_4(4) ] , KJ(1:half_KJ));

% LLR_EST_B0_right=log(pos_inv./neg_inv);
% LLR_EST_B0_center=log(Positive_1./negative_1);
% LLR_EST_B0_left=log(pos_inv_2./neg_inv_2);

% KJ_LEFT=KJ(1:half_KJ);
% KJ_RIGHT=KJ(half_KJ+1:end);

[~,c2]=min(abs(KJ-y));
ZR=LLR_FULL(c2);

% if ((y >= Threshold_left)&&(y <=Threshold_rigth))
%     [c1,c2]=min(abs(KJ-y));
%     ZR=LLR_EST_B0_center(c2)
% elseif( y<Threshold_left)
%     [d1,d2]=min(abs(KJ_LEFT-y));
%     ZR=LLR_EST_B0_left(d2)
% elseif ( y>Threshold_rigth)
%     u=2
%     [e1,e5]=min(abs(KJ_RIGHT-y));
%     ZR=LLR_EST_B0_right(e5)

end




