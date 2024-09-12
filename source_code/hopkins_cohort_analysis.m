close all; clc; clear
runtag = '240911';
%% PART 1: LOADING Morphometric Data
xlsname_datametrics='data.xlsx';
xlsname_featureInfo='feature_list.xlsx';
[num,~,all2]=xlsread(xlsname_datametrics,'by section all');
ccsi=strcmpi(all2(3,:),'SI'); %find columns with sample information
numlbl=all2(2:4,:); %column labels
% create data metrics
simat=num(:,ccsi); % sample info metrics; NaN for string 
simatc=all2(:,ccsi); %5th row and on contains data matrix
fmat=num(:,~ccsi); % feature metrics
silbl=numlbl(:,ccsi);
flbl=numlbl(:,~ccsi);
flbl=[num2cell(1:size(flbl,2));flbl]; % adding  identifyinig number in first row 

% sample info label; add one more column called sectionnum
silbl(3,end+1)={'SectionNum'};
silbl(2,end)={'SI'};
silbl(1,end)={nan}; 
silbl=[num2cell(1:size(silbl,2));silbl]; % adding identifyinig number
% ensemble patient data from sections.
cck=strcmpi(silbl(end,:),'patient'); %find patient ID column
up=unique(simat(:,cck));
up(isnan(up)) = [];
% up=up(isback&iscluerange);
plm=[]; pls=[]; pli=[]; plmed=[]; % patient level (PL) metrics, mean/std/info/median
plic={};
% iterate unique patient id
for kp=1:length(up)
    cc=simat(:,cck)==up(kp);
    pm=nanmean(fmat(cc,:),1);
    pmed=nanmedian(fmat(cc,:),1);
    ps=nanstd(fmat(cc,:),[],1);
    pi0=simat(cc,:);
    pi0=pi0(1,:);
    pi0c=simatc(cc,:);
    pi0c=pi0c(1,:);
    plm=[plm;pm]; % patient level mean (among sections)  
    plmed=[plmed;pmed]; % patient level median (among sections)
    pls=[pls;ps]; % patient level std (among secionts)
    pli=[pli;[pi0,sum(cc)]];
    plic=[plic;[pi0c,{sum(cc)}]];
end

%% OUTPUT 1: median value of morphometric features from tissue sections from each skin donor
save('plmed.mat',"plmed");
%% PART 2: Sample cohort selection
isMultiSec = cell2mat(plic(:,end))>1;

cc=1:size(flbl,2);
[cgender]=find(strcmpi(silbl(end,:),'gender'),1,'first');
ismale=strcmpi(plic(:,cgender),'male');
ccm=ismale; ccf=~ccm;

% select body part
[cbodypart]=find(strcmpi(silbl(end,:),'body part'),1,'first');
isback=strcmpi(plic(:,cbodypart),'back'); 

% select race
[crace]=find(strcmpi(silbl(end,:),'race'),1,'first');
iswhite=strcmpi(plic(:,crace),'white');

% white and back body part selected for further analysis
ccsample = iswhite & isback;

% age 
[cage]=find(strcmpi(silbl(end,:),'age')==1); %column id of age
age = pli(:,cage);

%% OUTPUT 2: patient information
save("ismale.mat","ismale");save("iswhite.mat","iswhite");save("isback.mat","isback");save("age.mat","age")
%% load tissue feature info 
[~,~,allf]=xlsread(xlsname_featureInfo);
allflbl=allf(1,:);
allf(1,:)=[];
cnonsense=strcmpi(allflbl,'nonsense');
isnonsense=cell2num(allf(:,cnonsense));
%% FIGURE 2C: Correlation with aging and Cohen's D
cc=pli(:,end)>1;
plcv=abs(pls(:,:)./plm(:,:)); 
plcv(~cc,:)=nan;

% setup criteron matrix and generate correlation matrix for each criteria
crimat{1}=plm>-inf; 
crimat{2}=plcv<0.5;  % CV below 0.5
crimat{3}=plcv<1;  % CV below 1

criname{1}='all';
criname{2}='CVless50';
criname{3}='CVless100';
    
lbltmp0={'N-all','N-male','N-female','rho-all', 'rho-male', 'rho-female',...
    'Pval-all','Pval-male','Pval-female','cohend-all','cohend-male','cohend-female'};

plmuse=plmed;

cc=1:size(plmuse,2);

restmp0=[];
restmp0lbl={};
restmp0lbl2={};
tic;
% Iterate criterion matrix
for kss=1:length(crimat)
    cef=0; cefm=0; ceff=0;
    p=0; pm=0; pf=0;
    Nm=0; Nf=0; Nall=0;
    cohendall=0;
    cohendm=0;
    cohendf=0;
    % Iterate each feature column
    for kf=1:length(cc)
          % exclude nan value in the colum 
           tst= ccsample;
           tst= tst & ~isnan(plmuse(:,cc(kf)));
           tst= tst & crimat{kss}(:,cc(kf)); % adding set -criterion
           % can add more criteron such as stainning batch.. or exlcude
           % outlier
           Nall(kf)=sum(tst);
           Nm(kf)=sum(tst& ccm);
           Nf(kf)=sum(tst& ccf);

            if sum(tst>0)
                ageuse = pli(tst,cage);
                plmusek = plmuse(tst,cc(kf));
               [cef(kf),p(kf)]=corr(ageuse,plmusek,'type','spearman');
                isold0 = ageuse>70; isyoung0 = ageuse<30;
                matA = plmusek(isold0,:);matB = plmusek(isyoung0,:);
                cohendall(kf) = cohend(matA,matB);
            else
                cef(kf)=0;
                p(kf)=0;
                cohendall(kf)=0;
            end

            if sum(tst & ccm>0)
                ageuse = pli(tst&ccm,cage);
                plmusek = plmuse(tst & ccm,cc(kf));
               [cefm(kf),pm(kf)]=corr(ageuse,plmusek,'type','spearman');
                isold0 = ageuse>70; isyoung0 = ageuse<30;
                matA = plmusek(isold0,:);matB = plmusek(isyoung0,:);
                cohendm(kf) = cohend(matA,matB);
            else
                cefm(kf)=0;
                pm(kf)=0;
                cohendm(kf)=0;
            end

           if sum(tst&ccf>0)
               ageuse = pli(tst&ccf,cage);
               plmusek = plmuse(tst & ccf,cc(kf));
               [ceff(kf),pf(kf)]=corr(ageuse,plmusek,'type','spearman');
                isold0 = ageuse>70; isyoung0 = ageuse<30;
                matA = plmusek(isold0,:);matB = plmusek(isyoung0,:);
                cohendf(kf) = cohend(matA,matB);
           else
               ceff(kf)=0;
               pf(kf)=0;
               cohendf(kf)=0;
           end            
    end
    restmp0=[restmp0, [Nall' Nm' Nf' cef' cefm'  ceff' p' pm' pf' cohendall' cohendm' cohendf']]; %ceff and pf are zero
    restmp0lbl=[restmp0lbl,lbltmp0];
    restmp0lbl2=[restmp0lbl2,repmat(criname(kss),1,length(lbltmp0))];
end
restmp=[flbl(:,cc)' num2cell(restmp0)];
reslbltmp=cell(1,size(restmp,2));
reslbltmp(1,size(flbl,1)+1:end)=restmp0lbl2;
reslbltmp(2,size(flbl,1)+1:end)=restmp0lbl; 
if 1
    xlsname_corrmet=['Feature Age correlation metrics',runtag,'.xlsx'];
    xlswrite(xlsname_corrmet,restmp,'correlation','A3');    
    xlswrite(xlsname_corrmet,reslbltmp,'correlation','A1');
end
toc;
%% OUTPUT 3: selected age associated features of white back based on FIGURE 2C
fall=cell2mat(restmp(:,[8 11 14])); %8:rho 11:p-val 14:cohend
age_associated_features=abs(fall(:,1))>0.3 & fall(:,3)>1 & ~isnonsense;
save('age_associated_features.mat',"age_associated_features");
%% FIGURE 2C: 
figure;plot(fall(~age_associated_features& ~isnonsense,1),fall(~age_associated_features& ~isnonsense,3),'k.');hold on;plot(fall(age_associated_features,1),fall(age_associated_features,3),'r.');bjff3;xlim([-1 1]);ylim([0 3]);xticks([-1:0.5:1]);yticks([0:0.5:3])
%% Figure 2F
% selected sample and features
pwset = plmuse(ccsample,:);
pwset(isnan(pwset))=0; %replace inf to zero 
ageuse = age(ccsample);
[agesorted,I]=sort(ageuse);
pwsetS = pwset(I,age_associated_features);
% pairwise correlation
[rho,pval] = corr(pwsetS,'Type','Spearman','Rows','pairwise');

cg2 = clustergram(rho,'Colormap',redbluecmap,'RowPdist', 'correlation',...
                             'ColumnPdist', 'correlation');
set(cg2,'Dendrogram',0.32) 
set(gcf,'Position',[100 100 413 400])
%% Write table with core process ID to describe core in plain english
%Feature ID of selected features
ogLUT= find(ccf); 

% define c1~c8 from clustergram manually here before proceeding to next one
cla = {c1,c2,c3,c4,c5,c6,c7,c8,c9,c10};
tmps = [];
% for each cluster, link feature ID to a core ID
for i=1:length(cla)
tmp = cellfun(@str2num, cla{i}.RowLabels); %convert sfID to number
tmpog = ogLUT(tmp); %query feature ID from sfID
tmp2 = [tmp,tmpog,ones([length(tmp),1]).*i]; %sfID, feature ID, cluster ID
tmps = [tmps;tmp2];
end
% sort based on sfID
[B,I]=sort(tmps(:,1));
tmpss = tmps(I,:);
% add feature info into the table
restmp=[flbl(2:end,tmpss(:,2))', num2cell(tmpss),  num2cell(fall(tmpss(:,2),:))];
restmp=[{'tissue component','feature group','variable','corridx','ogidx','clgroup','rho','p','d'};restmp];
xlsname_feat_cluster=['Feature_cluster',runtag,'.xlsx'];
writecell(restmp,xlsname_feat_cluster)


%%
fg = restmp(2:end,2);
fg = cellfun(@(x) str2double(x(2)), fg, 'UniformOutput', false);
fg = cell2mat(fg);
cg = restmp(2:end,6);
cg = cell2mat(cg);
fcg=[fg,cg];

%count occurence of feature group and core process pair
[C,ia,ic] = unique(fcg,'row');
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];

% Plot sankey diagram for Figure 2
figure;
% Sankey diagram panel height
data1{1}=[2 10 63 48]; 
data1{2}=[8 6 14 12 26 48 4 5];
data2{1}=data1{2};
% Panel width
w = 20; 
% Size of flow from first panels to second panels
data{1}=[0 0 0 0 2 0 0 0;8 0 0 0 0 2 0 0;0 1 0 5 24 25 4 4;0 5 14 7 0 21 0 1];
% x-axis scale
X=[0 10];
% panel color
cmap = [220 104 238;
    255 125 70;
    187 255 119; 
    255 199 32; 
    35 165 255;
    255 59 59;
    255 37 200;
    112 48 160];
cmap = cmap./255;

barcolors{1}=[244 67 54;156 39 176;205 220 57;255 84 34]./255;
barcolors{2}=cmap;
% flow color
c = [.7 .6 .3];
for j=1:1
    if j>1
        ymax=max(ymax,sankey_yheight(data1{j-1},data2{j-1}));
        y1_category_points=sankey_alluvialflow(data1{j}, data2{j}, data{j}, X(j), X(j+1), y1_category_points,ymax,barcolors{j},barcolors{j+1},w,c);
    else
        y1_category_points=[];
        ymax=sankey_yheight(data1{j},data2{j});
        y1_category_points=sankey_alluvialflow(data1{j}, data2{j}, data{j}, X(j), X(j+1), y1_category_points,ymax,barcolors{j},barcolors{j+1},w,c);
    end
end  
%% univariate raw
runtag = '240108_univariate';
plmusek = plmuse;
plmusek(isnan(plmusek))=0;
plmusek = plmusek(ccsample,:);
cfmat=plmusek; 
% cfmat=zscore(plmusek); % there is no effect doing zscore for univariate
y0=pli(ccsample,cage); 
ccm=ccm(ccsample);
pmresuf = univariateglm(cfmat,y0,runtag,flbl,ccm);
%mean(err) std(err) male-err female-err
%% FIG 3A: plot rho vs univaraite MAE of correlated features
MAEint = pmresuf(ccf,1);
rho = cef(ccf);
p = plot(rho,MAEint,'k.');bjff3;
p.MarkerSize = 10;
xlim([-0.9 0.9]);ylim([12 18])
%% FIG 3A: plot histogram
[tmp,edges]=histcounts(MAEint,[12:.1:18]);
bar(tmp);bjff3;ylim([0 15])
set(gca, 'XDir','reverse');xticks([]);set(gca,'YAxisLocation','right')
set(gcf, 'Position',[100 100 300 100])
%% Table 1A: univaraite error
restmp=[flbl(2:end,:)', num2cell([ccf,[1:1090]',pmresuf])];
restmp=[{'tissue component','feature group','variable','selected-feat','ogidx','mean-err','std-err','male-err','female-err'};restmp];
writecell(restmp,xlsname_feat_cluster,'Sheet',3)
%% FIG 3C: distribution along cores
unierr = MAEint;
ec1=unierr(clgroup==1);
ec2=unierr(clgroup==2);
ec3=unierr(clgroup==3);
ec4=unierr(clgroup==4);
ec5=unierr(clgroup==5);
ec6=unierr(clgroup==6);
ec7=unierr(clgroup==7);
ec8=unierr(clgroup==8);
ec9=unierr(clgroup==9);
ec10=unierr(clgroup==10);
disp('mean :')
disp([mean(ec1),mean(ec2),mean(ec3),mean(ec4),mean(ec5),mean(ec6),mean(ec7),mean(ec8),mean(ec9),mean(ec10)])
violin_dots({ec1,ec2,ec3,ec4,ec5,ec6,ec7,ec8,ec9,ec10});bjff3;set(gcf,'position',[100 100 720 300]);xticks([1:10]);xlim([0,11]);xticklabels([1:10]);ylim([11,20])
bjff3;xlim([0,11]);ylim([10,20]);xticks([1:10]);yticks([10,13,16,19]);
set(gcf,'Position',[100 100 350 350])

%% pair-wise ttest
ecs={ec1,ec2,ec3,ec4,ec5,ec6,ec7,ec8,ec9,ec10};
hs=zeros([10,10]);
ps=zeros([10,10]);
for i=1:10
    for j=1:10
        [h,p]=ttest2(ecs{i},ecs{j});
        hs(i,j)=h;
        ps(i,j)=p;
    end
end

%% FIG 3C: circular tree
treeList={[],[]};
k=1;
classNameSet={'A','B','C','D','E','F','G','H','I','J'};
for i=1:10
    for j=1:sum(clgroup==i)
        treeList{k,1}=[classNameSet{i}];
        treeList{k,2}=[num2str(sum(clgroup==i))];
        k=k+1;
    end
end
CT=circleTree2(treeList);
CT=CT.draw();
%manually delete all labels for now because I don't know how else

%% define top univariate features
% univaraite_err = [[1:1090]',pmresuf];
% univaraite_err_sorted = sortrows(univaraite_err,2);
% topN = 10; % how many number of features do you want?
% topN_uni_idx = univaraite_err_sorted(1:topN,1);
% topN_uni_err = sortrows(univaraite_err_sorted(1:topN,1:2),1);
% [sharedvals,idx] = intersect(tmpss(:,2),topN_uni_idx,'stable');
% core = [tmpss(idx,2:3),topN_uni_err]; %ogid,clid,mae
% core = sortrows(core,4);
% ccfeat_top_uni = zeros([1090,1]);
% ccfeat_top_uni(topN_uni_idx)=1;
% ccfeat_top_uni=logical(ccfeat_top_uni);
%% define top bivariate features
% ccfeat_top_bi = zeros([1090,1]);
% ccfeat_top_bi(a)=1;
% ccfeat_top_bi=logical(ccfeat_top_bi);
%% bivariate of correlated features
runtag = '240108';
ccfeat=ccf;
vnew = 1:sum(ccf);
C = nchoosek(vnew,2);

pmresuf=[];
cfmat=zscore(plmuse(ccsample,logical(ccfeat)));
%no PCA, we need specific pair of variables
y0=pli(ccsample,cage); 
for kff=1:length(C) %iterate each combination for bivariate 
    corfeatids = C(kff,:);
    x0=cfmat(:,corfeatids);
    yerrs=[];randyerrs=[];
    yts=[];ypreds=[];
    parfor kcv=1:length(y0) %iterate leave one out for samples
        p=[1:kcv-1,kcv+1:length(y0),kcv];
        trainnum=length(p)-1;
        trainid=p(1:trainnum);
        testid=p(trainnum+1:end);
        x=x0(trainid,:);
        y=y0(trainid);
        randy=y(randperm(length(y)));
        xt=x0(testid,:);
        yt=y0(testid);
        % glm
        mdl = fitglm(x,y);
        ypred = predict(mdl,xt);
        yerr=abs(ypred-yt);
        yerrs =[yerrs;yerr];
        % permutated glm
        randmdl = fitglm(x,randy);
        randypred = predict(randmdl,xt);
        randyerr=abs(randypred-yt);
        randyerrs = [randyerrs;randyerr];
        yts=[yts;yt];ypreds=[ypreds;ypred];
    end 
    ogfeatids = ogLUT(corfeatids)';
    clgroups = [tmpss(tmpss(:,2)==ogfeatids(1),3) tmpss(tmpss(:,2)==ogfeatids(2),3)];
    if length(clgroups)<2;clgroups=[clgroups,0];end
    pmresuf=[pmresuf;[ogfeatids clgroups mean(yerrs) mean(randyerrs) std(yerrs)]];
end
% note: find "bad" patient 
xlsname_bipred=['Bivariate Predicting power',runtag,'.xlsx'];
pmresuf = sortrows(pmresuf,5); %sort by mae

bipred=[num2cell(1:length(pmresuf(:,1)))',flbl(2:end,pmresuf(:,1))', flbl(2:end,pmresuf(:,2))', num2cell(pmresuf)];
bipred=[{'rank','tissue component A','feature group A','variable name A',...
    'tissue component B','feature group B','variable name B',...
    'variable id A','variable id B','clgroup A','clgroup B','MAE','permMAE','stdev_err'};bipred];
writecell(bipred,xlsname_bipred);
%% FIG 3F: plot bivariate pairs by clgroup
clpair = pmresuf(:,3:4);
%count occurence of feature group and core process pair
[C,ia,ic] = unique(clpair,'row');
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];

value_counts(value_counts(:,2)==0,:)=[];
Data = zeros(max(value_counts(:,1)),max(value_counts(:,2)));
for i = 1:size(value_counts,1)
    Data(value_counts(i,1),value_counts(i,2))=value_counts(i,3);
end
%Data = 8 x8 matrix  
myColorMap = zeros(length(Data),3);
% myColorMap = lines(length(Data));
circularGraph(Data,'Colormap',myColorMap);
%% FIG3D: Compare univaraite MAE and bivariate MAE
pmresuf = xlsread("Bivariate Predicting power240108_t1_bivariate_raw.xlsx");
univaraite_err = xlsread("Univariate Predicting power240108_univariate.xlsx");
univaraite_err = univaraite_err(:,[1,5]);
univaraite_err_sorted = sortrows(univaraite_err,2);

bimae = pmresuf(:,5);
meanunimae = [];
for i = 1:size(bimae,1)
    meanunimae=[meanunimae;mean([univaraite_err_sorted(univaraite_err_sorted(:,1)==pmresuf(i,1),2),univaraite_err_sorted(univaraite_err_sorted(:,1)==pmresuf(i,2),2)])];
end
plot(bimae,meanunimae,'k.');bjff3;axis equal;
xlim([11.5 17.5]);ylim([11.5 17.5])
hold on; plot([11:18],[11:18],'k--')

figure;
scale='log';
gradenum=100;
gradelimit=100;              
DensitySCPlotW(bimae,meanunimae,gradenum,gradelimit,scale);bjff3;axis equal;
xlim([11.5 17.5]);ylim([11.5 17.5])
hold on; plot([11:18],[11:18],'k--')
xlabel('Bivariate MAE');ylabel('Mean univariate MAE of the pair');
%% bivariate linear regression

% manually define bivariate pair of interest ex)4 and 27 here
x1=plmuse(:,4);x2=plmuse(:,27);y=ageuse;

X = [ones(size(x1)) x1 x2 x1.*x2];
b = regress(y,X);
scatter3(x1,x2,y,'filled')
hold on
x1fit = min(x1):0.01:max(x1);
x2fit = min(x2):0.01:max(x2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
mesh(X1FIT,X2FIT,YFIT); zlim([0,100]); bjff3;
xlabel('Hair composition [%]')
ylabel('Anisotropy_B_masked [-]')
zlabel('Age [years]')
%% tri-variate
% runtag = '230117_t1_trivariate_raw';
% ccfeat=ccfeat_global;
% vnew = 1:sum(ccfeat_global);
% C = nchoosek(vnew,3);
% 
% pmresuf=[];
% cfmat=zscore(plmuse(ccsample,logical(ccfeat)));
% %no PCA, we need specific pair of variables
% y0=pli(ccsample,cage); 
% for kff=1:length(C) %iterate each combination for bivariate 
%     if mod(kff,100) == 0;fprintf('At iteration %d...\n',kff);end
%     corfeatids = C(kff,:);
%     x0=cfmat(:,corfeatids);
%     yerrs=[];randyerrs=[];
%     parfor kcv=1:length(y0) %iterate leave one out for samples
%         p=[1:kcv-1,kcv+1:length(y0),kcv];
%         trainnum=length(p)-1;
%         trainid=p(1:trainnum);
%         testid=p(trainnum+1:end);
%         x=x0(trainid,:);
%         y=y0(trainid);
%         randy=y(randperm(length(y)));
%         xt=x0(testid,:);
%         yt=y0(testid);
%         % glm
%         mdl = fitglm(x,y);
%         ypred = predict(mdl,xt);
%         yerr=abs(ypred-yt);
%         yerrs =[yerrs;yerr];
%         % permutated glm
%         randmdl = fitglm(x,randy);
%         randypred = predict(randmdl,xt);
%         randyerr=abs(randypred-yt);
%         randyerrs = [randyerrs;randyerr];
%     end 
%     ogfeatids = ogLUT(corfeatids)';
% %     clgroups = [tmpss(tmpss(:,2)==ogfeatids(1),3) tmpss(tmpss(:,2)==ogfeatids(2),3)];
% %     if length(clgroups)<2;clgroups=[clgroups,0];end
%     pmresuf=[pmresuf;[ogfeatids mean(yerrs) mean(randyerrs) std(yerrs)]];
% end
% % note: find "bad" patient 
% xlsname_tripred=['Trivariate Predicting power',runtag,'.xlsx'];
% pmresuf = sortrows(pmresuf,4); %sort by mae
% 
% tripred=[num2cell(1:length(pmresuf(:,1)))',flbl(2:end,pmresuf(:,1))', flbl(2:end,pmresuf(:,2))',flbl(2:end,pmresuf(:,3))', num2cell(pmresuf)];
% tripred=[{'rank','tissue component A','feature group A','variable name A',...
%     'tissue component B','feature group B','variable name B',...
%     'tissue component C','feature group C','variable name C',...
%     'variable id A','variable id B','variable id C','MAE','permMAE','stdev_err'};bipred];
% writecell(tripred,xlsname_tripred);
%% preprocess input

% Input: plmed, ccsample, ccfeat,age
load('plmed.mat')
load('isback.mat')
load('iswhite.mat')
ccsample = isback&iswhite;
load('ccf.mat')
load('age.mat')
y0=age(ccsample); 
% Output: refined plmed, age
runtag = '230125_Npc'; 
ccsample= isback&iswhite;
% ccfeat_use = ccfeat_top_uni;
ccfeat_use = ccf;
plmusek = plmed;
plmusek(isnan(plmusek))=0;
plmusek = plmusek(ccsample,ccfeat_use);
cfmat=zscore(plmusek);
dopca=1;
if dopca
[coeff,score,latent]=pca(cfmat);
% save('hopkins_old_ccfeatglobal_PCAcoeff.mat','coeff');
vcum=cumsum(latent)/sum(latent);
v95num=find(vcum-0.95>0,1,'first'); %33 for white back; 34 for multisec
disp(sum(v95num)) %display 
x0=score(:,1:v95num);  
else
x0=cfmat;
end
% xg0=[x0 ismale(ccsample)]; % including gender for building model
% ismalignant2 = ismalignant(ccsample);
%% FIG3H: PC1 vs PC2 for age groups
% ageuse=y0;
ageyoung = ageuse<30;
agemiddle = 30<=ageuse & ageuse<=60;
ageold = 60<ageuse;
[sum(ageyoung),sum(agemiddle),sum(ageold)]
figure;plot(x0(ageyoung,1),x0(ageyoung,2),'r.');hold on
plot(x0(agemiddle,1),x0(agemiddle,2),'g.');
plot(x0(ageold,1),x0(ageold,2),'b.');bjff3;
%% number of input variable VS. accuracy
maess=[];errmats=[];
x00=score;
for i=1:size(x00,2)
    x0=x00(:,1:i);
    errmat=[];ypredall=[];
    parfor kcv=1:length(y0)
        p=[1:kcv-1,kcv+1:length(y0),kcv];
        trainnum=length(p)-1;
        trainid=p(1:trainnum);
        testid=p(trainnum+1:end);
        x=x0(trainid,:);
        y=y0(trainid);
        randy = y(randperm(length(y)));
        xt=x0(testid,:);
        yt=y0(testid);
        % glm
        mdl = fitglm(x,y);
        ypred = predict(mdl,xt);
        yerr=(ypred-yt);
        randmdl = fitglm(x,randy);
        randypred = predict(randmdl,xt);
        randyerr=(randypred-yt);
        % lasso
        lambda = 1e-03;
        [B, FitInfo] = lasso(x,y,'Alpha',0.5,'CV',10);
        idxLambda1SE = FitInfo.Index1SE;
        idxLambda1SE = FitInfo.IndexMinMSE;
        coef = B(:,idxLambda1SE);
        coef0 = FitInfo.Intercept(idxLambda1SE);
        ypred2 = xt*coef + coef0;
        yerr2= (ypred2-yt);
        [B, FitInfo] = lasso(x,randy,'Alpha',0.5,'CV',10);
        idxLambda1SE = FitInfo.Index1SE;
        idxLambda1SE = FitInfo.IndexMinMSE;
        coef = B(:,idxLambda1SE);
        coef0 = FitInfo.Intercept(idxLambda1SE);
        randypred2 = xt*coef + coef0;
        randyerr2= (randypred2-yt);
        % svm
        svmmdl = fitrsvm(x,y,'KernelFunction','linear','Epsilon',0.1,'Alpha',ones(size(y,1),1)*0.001 );
        ypred3 = predict(svmmdl,xt);
        yerr3=(ypred3-yt);
        randsvmmdl = fitrsvm(x,randy,'KernelFunction','linear','Epsilon',0.1,'Alpha',ones(size(y,1),1)*0.001 );
        randypred3 = predict(randsvmmdl,xt);
        randyerr3=(randypred3-yt);
        %kernel
        kernelmdl = fitrkernel(x,y,'Learner','svm','NumExpansionDimensions',2048,'KernelScale',23.8858,'Lambda',2.9317e-04,'Epsilon',2.3795);
        randMdl = fitrkernel(x,randy,'Learner','svm','NumExpansionDimensions',2048,'KernelScale',23.8858,'Lambda',2.9317e-04,'Epsilon',2.3795);
        ypred4 = predict(kernelmdl,xt);
        randypred4 = predict(randMdl,xt);
        yerr4=mean(abs(ypred4-yt));
        randyerr4=mean(abs(randypred4-yt));
        errmat=[errmat;[testid(:) yerr randyerr yerr2 randyerr2 yerr3 randyerr3 yerr4 randyerr4]];
    end
    maes = mean(abs(errmat(:,2:end)));
    randcorr =[corr(errmat(:,2),errmat(:,3),'type','pearson'),corr(errmat(:,3),errmat(:,4),'type','pearson'),corr(errmat(:,6),errmat(:,7),'type','pearson'),corr(errmat(:,8),errmat(:,9),'type','pearson')];
    errmats=[errmats;errmat];
    maess = [maess;[maes,randcorr]];
end
maes = [[1:size(x00,2)]',maess(:,1:2),maess(:,9),maess(:,3:4),maess(:,10),maess(:,5:6),maess(:,11),maess(:,7:8),maess(:,12)];
% maes = [maes;mean(maes)];
maes_header = {'iterations','GLM','GLM_perm','GLM_corr','Lasso','Lasso_perm','Lasso_corr','SVM','SVM_perm','SVM_corr','Kernel','Kernel_perm','Kernel_corr'};
maes2 = [maes_header;num2cell(round(maes,2))];
xlsname_permute=['Npc_vs_MAE_',runtag,'.xlsx'];
writecell(maes2,xlsname_permute);

plot(maes(:,1),maes(:,2),'-');hold on; plot(maes(:,1),maes(:,5),'-'); plot(maes(:,1),maes(:,8),'-');plot(maes(:,1),maes(:,11),'-');bjff3
legend('GLM','Lasso','SVM','Kernel')
%% FIG 3I multivariate regression
% select Ndarray = 29D for here
root = '\\10.99.68.51\kyu\skin_aging\source_code\streamline\modeling\230117';
path(path,root);
load(fullfile(root,'plmed.mat'));
load(fullfile(root,'isback.mat'));
load(fullfile(root,'iswhite.mat'));
load(fullfile(root,'ccf.mat'));
load(fullfile(root,'age.mat'));
load(fullfile(root,'ismale.mat'));


ccsample = isback&iswhite;
ccm=ismale(ccsample);
y0=age(ccsample);

plmusek = plmed;
plmusek(isnan(plmusek))=0;
plmusek = plmusek(ccsample,ccf);
cfmat=zscore(plmusek);
[coeff,score,latent]=pca(cfmat);

vcum=cumsum(latent)/sum(latent);
v95num=find(vcum-0.95>0,1,'first'); %33 for white back; 34 for multisec
disp(sum(v95num)) %display 

%plot 95 variance
% figure;plot(0:1:10,[0,vcum],'o-');bjff3;xlim([0 11]);ylim([0 1.1])

x0=score(:,1:v95num);  
m=[x0,y0];
writematrix(m,'\\10.99.68.51\Kyu\skin_aging\source_code\hopkins_regression\M.csv');

rtree = fitrtree(x0,y0,'Leaveout','on'); %tunelength5, 100tree,
ypredr = predict(rtree,x0(1,:));
yerrr=(ypredr-y0);


[maes,errmat,y_glm,y_las,y_svm,y_rtree,randcorr] = agingmodel(x0,y0); %glm, lasso, svm
% check if randcorr is low enough, otherwise the model is questionable
disp(maes);

% output trained SVM for external dataset
svmmdl = fitrsvm(x0,y0,'KernelFunction','linear','Epsilon',0.1,'Alpha',ones(size(y0,1),1)*0.001 );
svmpred = predict(svmmdl,x0);
svmerr=(svmpred-y0);  %MAE=5.8155 with all dataset
disp(mean(abs(svmerr)));
save("svmmdl.mat",'svmmdl');

aes = abs(errmat(:,2:end));
figure;plot(y0,y_svm,'k.');bjff3;ylim([0 100]);xticks([0:20:100]);yticks([0:20:100]);hold on;plot(1:100,1:100,'k--');xticks([0:20:100]);xticklabels([])
figure;violin_dots({aes(:,1),aes(:,3),aes(:,5),tmp});bjff3;xticks([]);ylim([0,40]);yticks([0:10:40]);xticks([0:20:100]);xticklabels([]);

figure;plot(ry(rccm),rpred(rccm),'b+');hold on;plot(ry(~rccm),rpred(~rccm),'ro');plot(1:100,1:100,'k--');ylim([0 100]);xticks([0:20:100]);yticks([0:20:100]);xticks([0:20:100]);xticklabels([]);yticklabels([]);bjff3

% ttest to see SVM has significantly lower MAE
[h,p]=ttest2(aes(:,1),aes(:,3));disp([h,p])
[h,p]=ttest2(aes(:,1),aes(:,5));disp([h,p])
[h,p]=ttest2(aes(:,3),aes(:,5));disp([h,p])
%% FIG 4A: multivariate MAE by gender
figure;plot(y0(ccm),y_svm(ccm),'b+');hold on;plot(y0(~ccm),y_svm(~ccm),'ro');bjff3;ylim([0 100]);xticks([0:20:100]);yticks([0:20:100]);plot(1:100,1:100,'k--');xticklabels([]);
figure;violin_dots({abs(errmat(ccm,6)),abs(errmat(~ccm,6))});bjff3;ylim([0,25]);xticks([])
%% FIG 4B: aging feature labeled by gender
gen_labels = [];
genstat = [];
for i=1:108
    x_obs = cfmat(:,i);
    [cm,pm] = corr(x_obs(ccm2),y0(ccm2)); 
    [cf,pf] = corr(x_obs(~ccm2),y0(~ccm2));
    genstat=[genstat;[cm,pm,cf,pf]];
    if 0.3<=abs(cm) && pm<=0.05 && 0.3<=abs(cf) && pf<=0.05
        gen_label=3;
    elseif 0.3<=abs(cf) && pf<=0.05
        gen_label=2;
    elseif 0.3<=abs(cm) && pm<=0.05
        gen_label=1;
    else 
        gen_label=0;
    end
    gen_labels = [gen_labels,gen_label];
end
disp([sum(gen_labels==0),sum(gen_labels==1),sum(gen_labels==2),sum(gen_labels==3)])
genstat=[genstat,gen_labels'];
%% FIG 4C change in predicted age with gender in the input feature

%OPTION 1: add sex before pca (For now,  I like this)
plmuseg = [plmusek,ccm]; %add sex
cfmat=zscore(plmuseg);
[coeff,score,latent]=pca(cfmat);

vcum=cumsum(latent)/sum(latent);
v95num=find(vcum-0.95>0,1,'first'); %33 for white back; 34 for multisec
disp(sum(v95num)) %display 
x0=score(:,1:v95num);  %pc=29 without sex, pc=30 with sex
[maes1,errmat1,y_glm,y_las,y_svm1,randcorr] = agingmodel(x0,y0); %glm, lasso, svm
% check if randcorr is low enough, otherwise the model is questionable
% disp(randcorr);
disp(maes1)
aes1 = abs(errmat1(:,2:end));
figure;plot(y0(ccm),y_svm1(ccm),'b+');hold on;plot(y0(~ccm),y_svm1(~ccm),'ro');bjff3;ylim([0 100]);xticks([0:20:100]);yticks([0:20:100]);plot(1:100,1:100,'k--');xticklabels([]);yticklabels([]);
% figure;violin_dots({aes1(:,1),aes1(:,3),aes1(:,5)});bjff3;xticks([]);ylim([0,40]);yticks([0:10:40]); 
ed1=aes1(:,5)-aes(:,5); %change in absolute error
figure;plot(y0(ccm),ed1(ccm),'b+');hold on;plot(y0(~ccm),ed1(~ccm),'ro');bjff3;xlim([0 100]);plot(1:100,zeros(100,1),'k--');
disp([mean(ed1),mean(ed1(ccm)),mean(ed1(~ccm))])
ed2=errmat1(:,6)-errmat(:,6); %change in error
figure;plot(y0(ccm),ed2(ccm),'b+');hold on;plot(y0(~ccm),ed2(~ccm),'ro');bjff3;xlim([0 100]);plot(1:100,zeros(100,1),'k--');
disp([mean(ed2),mean(ed2(ccm)),mean(ed2(~ccm))])
%OPTION 2: add sex after pca and zscore 
cfmat=zscore(plmusek);
[coeff,score,latent]=pca(cfmat);

vcum=cumsum(latent)/sum(latent);
v95num=find(vcum-0.95>0,1,'first'); %33 for white back; 34 for multisec
disp(sum(v95num)) %display 
x0=score(:,1:v95num);  %pc=29 without sex, pc=30 with sex
x0=[x0,ccm]; %add sex
x0=zscore(x0);

[maes2,errmat2,y_glm,y_las,y_svm2,randcorr] = agingmodel(x0,y0); %glm, lasso, svm
% check if randcorr is low enough, otherwise the model is questionable
% disp(randcorr);
disp(maes2)
aes2 = abs(errmat2(:,2:end));
% figure;plot(y0(ccm),y_svm2(ccm),'b+');hold on;plot(y0(~ccm),y_svm2(~ccm),'ro');bjff3;ylim([0 100]);xticks([0:20:100]);yticks([0:20:100]);plot(1:100,1:100,'k--');
% figure;violin_dots({aes2(:,1),aes2(:,3),aes2(:,5)});bjff3;xticks([]);ylim([0,40]);yticks([0:10:40]); 
ed3=aes2(:,5)-aes(:,5);
figure;plot(y0(ccm),ed3(ccm),'b+');hold on;plot(y0(~ccm),ed3(~ccm),'ro');bjff3;xlim([0 100]);plot(1:100,zeros(100,1),'k--');
disp([mean(ed3),mean(ed3(ccm)),mean(ed3(~ccm))])
ed4=errmat2(:,6)-errmat(:,6); %change in error
figure;plot(y0(ccm),ed4(ccm),'b+');hold on;plot(y0(~ccm),ed4(~ccm),'ro');bjff3;xlim([0 100]);plot(1:100,zeros(100,1),'k--');
disp([mean(ed4),mean(ed4(ccm)),mean(ed4(~ccm))])
%% FIG 4D leave-one-out gender classification by N features


plmusek = plmed;
plmusek(isnan(plmusek))=0;
plmusek = plmusek(ccsample,ccf);
cfmat=zscore(plmusek);
[coeff,pcscore,latent]=pca(cfmat);
vcum=cumsum(latent)/sum(latent);
v95num=find(vcum-0.95>0,1,'first'); %33 for white back; 34 for multisec
disp(sum(v95num)) %display 
x0=pcscore(:,1:v95num);
Nfeat=size(x0,2); 
% x0=cfmat;
% Nfeat=29;
y0 = ccm;
I = rankfeatures(x0',y0,NumberOfIndices=Nfeat);
crs=[];rrs=[];
for i = 1:Nfeat
    yerrs=[];randyerrs=[];
    x00=x0(:,I(1:i));
% leave one out validation
parfor kcv=1:length(y0)
    p=[1:kcv-1,kcv+1:length(y0),kcv];
    trainnum=length(p)-1;
    trainid=p(1:trainnum);
    testid=p(trainnum+1:end);   
    x=x00(trainid,:);
    y=y0(trainid);
    randy = y(randperm(length(y)));
    xt=x00(testid,:);
    yt=y0(testid);
    % Discriminant analysis classifer
    ypred = classify(xt,x,y);
    yerr=abs(ypred-yt);
    randypred = classify(xt,x,randy);
    randyerr=abs(randypred-yt);
    randyerrs=[randyerrs;randyerr];yerrs=[yerrs;yerr];
end 
cr = 1-sum(yerrs)/length(y0);
rr = 1-sum(randyerrs)/length(y0);
crs=[crs;cr];
rrs=[rrs;rr];
end
disp(crs)
figure;plot(1:Nfeat,crs,'ko-');bjff3;
%% FIG 4D leave-one-out gender classification by N features - by accuracy
plmusek = plmed;
plmusek(isnan(plmusek))=0;
plmusek = plmusek(ccsample,ccf);
cfmat=zscore(plmusek);
x0=cfmat;
Nfeat=size(x0,2); 
% x0=cfmat;
% Nfeat=29;
y0 = ccm;
crs=[];rrs=[];
Nfeat = 10; %Nfeat cannot be more than observations -> Nsample=99
for i = 1:Nfeat
    yerrs=[];randyerrs=[];
    x00=x0(:,I(1:i));
% leave one out validation
parfor kcv=1:length(y0)
    p=[1:kcv-1,kcv+1:length(y0),kcv];
    trainnum=length(p)-1;
    trainid=p(1:trainnum);
    testid=p(trainnum+1:end);   
    x=x00(trainid,:);
    y=y0(trainid);
    randy = y(randperm(length(y)));
    xt=x00(testid,:);
    yt=y0(testid);
    % Discriminant analysis classifer
    ypred = classify(xt,x,y);
    yerr=abs(ypred-yt);

    % Generalized Linear Model Classifier (GLM)
    b = glmfit(x, y, 'binomial', 'link', 'logit','weights',w0);
    ypred = glmval(b, xt, 'logit');
    yerr = round(abs(ypred-yt));
    
    % GLM - random model
    br = glmfit(x, randy, 'binomial', 'link', 'logit','weights',w0);
    randypred = glmval(br, xt, 'logit');
    randyerr=round(abs(randypred-yt));

    % Linear discriminant classifier (LDA)
    % ypred = classify(xt,x,y);
    % yerr = abs(ypred-yt);
    % randypred = classify(xt,x,randy);
    % randyerr=abs(randypred-yt);

    randyerrs=[randyerrs;randyerr];yerrs=[yerrs;yerr];
end 
cr = 1-sum(yerrs)/length(y0);
rr = 1-sum(randyerrs)/length(y0);
crs=[crs;cr];
rrs=[rrs;rr];
end
disp(crs)
figure;plot(1:Nfeat,crs*100,'ko-');bjff3;
%%  FIG 4F: PC vs accuracy
top10=[279
836
697
543
37
296
302
30
36
330
];
plmusek = plmed;
plmusek(isnan(plmusek))=0;
plmusek = plmusek(ccsample,ccf);
cfmat=zscore(plmusek);
[coeff,score,latent]=pca(cfmat);
vcum=cumsum(latent)/sum(latent);
v95num=find(vcum-0.95>0,1,'first'); %33 for white back; 34 for multisec
disp(sum(v95num)) %display 
x0=score(:,1:v95num);   
y0 = ccm;
W = ones(1,99);
crs=[];rrs=[];
for i = 1:size(x0,2)
    yerrs=[];randyerrs=[];
    % x00=x0(:,1:i);
    x00=x0(:,i);

% leave one out validation
parfor kcv=1:length(y0)
    p=[1:kcv-1,kcv+1:length(y0),kcv];
    trainnum=length(p)-1;
    trainid=p(1:trainnum);
    testid=p(trainnum+1:end);   
    x=x00(trainid,:);
    y=y0(trainid);
    w0=W(trainid);
    randy = y(randperm(length(y)));
    xt=x00(testid,:);
    yt=y0(testid);
    % Discriminant analysis classifer
    % ypred = classify(xt,x,y);
    % yerr=abs(ypred-yt);
    % randypred = classify(xt,x,randy);
    % randyerr=abs(randypred-yt);
    % 
    % 
    % Generalized Linear Model Classifier
    b = glmfit(x, y, 'binomial', 'link', 'logit','weights',w0);
    ypred = glmval(b, xt, 'logit');
    yerr = round(abs(ypred-yt));
    br = glmfit(x, randy, 'binomial', 'link', 'logit','weights',w0);
    randypred = glmval(br, xt, 'logit');
    randyerr=round(abs(randypred-yt));
    
    randyerrs=[randyerrs;randyerr];yerrs=[yerrs;yerr];
end 
cr = 1-sum(yerrs)/length(y0);
rr = 1-sum(randyerrs)/length(y0);
crs=[crs;cr];
rrs=[rrs;rr];
end
disp(crs)
figure;plot(1:size(x0,2),crs*100,'ko-');bjff3;
%% FIG 4D: add sex to classify age again and show it does not improve
y0=age(ccsample);
plmusek = plmed;
plmusek(isnan(plmusek))=0;
plmusek = plmusek(ccsample,ccf);
cfmat=zscore(plmusek);
[coeff,score,latent]=pca(cfmat);

vcum=cumsum(latent)/sum(latent);
v95num=find(vcum-0.95>0,1,'first'); %33 for white back; 34 for multisec
disp(sum(v95num)) %display 


%define ccm value
ccm=ismale(ccsample);
% ccm=ccm-.5;
ccm=(ccm-0.5)*2*0;

x0=score(:,1:v95num);
x0=[x0,ccm]; %add sex as a feature

[maes,errmat2,y_glm2,y_las2,y_svm2,randcorr2] = agingmodel(x0,y0);
disp(maes)

% 7.95 [-50 50]
% 7.97 [-5 5]
% 7.99 [-2 2]
% 7.92 [-1 1]
% 7.86 [-.5 5] single
% 7.84 [-.5 5] double
% 7.74 [-.1 1]
% 7.72 [0 0]


% check if randcorr is low enough, otherwise the model is questionable
disp(randcorr);

aes2 = abs(errmat2(:,2:end));
figure;plot(y0,y_svm2,'k.');bjff3;ylim([0 100]);xticks([0:20:100]);yticks([0:20:100]);hold on;plot(1:100,1:100,'k--')
figure;violin_dots({aes2(:,1),aes2(:,3),aes2(:,5)});bjff3;xticks([]);ylim([0,40]);yticks([0:10:40]);

svm_d = y_svm2-y_svm;
figure;plot(y0(ccm2),svm_d(ccm2),'b+');hold on;plot(y0(~ccm2),svm_d(~ccm2),'r+');bjff3;ylim([-1.5 1.5])
%% FIG 4E: best univariate sex discriminant
xx=plmusek(:,I(1));
figure;plot(y0(ccm2),xx(ccm2),'b.');hold on; plot(y0(~ccm2),xx(~ccm2),'r.');bjff3;ylim([-30 360]);xticks([0:20:100]);

ym=xx(ccm2&ageyoung);
yf=xx(~ccm2&ageyoung);
om=xx(ccm2&ageold);
of=xx(~ccm2&ageold);
figure;violin_dots({ym,yf,om,of});bjff3;ylim([-30 360]);xticks([]);yticks([]);
%% FIG 4G: top 10 sex discriminant

load('plmed.mat')
load('isback.mat')
load('iswhite.mat')
load('ccf.mat')
load('age.mat')
load('ismale.mat')
ccsample = iswhite & isback;
ccm = ismale(ccsample);

plmusek = plmed;
plmusek(isnan(plmusek))=0;
plmusek = plmusek(ccsample,ccf);
cfmat=zscore(plmusek);

% option 1: by p-value
[I,Z] = rankfeatures(cfmat',ccm,NumberOfIndices=97);
ogLUT = find(ccf);
sfeatID = ogLUT(I);
% top10sd = flbl(2:end,sfeatID)';

% option 2: by LDA accuracy


fid0=279;
x=age(ccsample);
y=plmed(ccsample,fid0);
plot(x(ccm),y(ccm),'b+');hold on;plot(x(~ccm),y(~ccm),'ro');

%%

xlsname_datametrics='overall_220715.xlsx';
xlsname_featureInfo='feature_list_231107.xlsx';
[num,~,all2]=xlsread(xlsname_datametrics,'by section all');
ccsi=strcmpi(all2(3,:),'SI'); %find columns with sample information
numlbl=all2(2:4,:); %column labels
% create data metrics
simat=num(:,ccsi); % sample info metrics; NaN for string 
simatc=all2(:,ccsi); %5th row and on contains data matrix
fmat=num(:,~ccsi); % feature metrics
silbl=numlbl(:,ccsi);
flbl=numlbl(:,~ccsi);
flbl=[num2cell(1:size(flbl,2));flbl]; % adding  identifyinig number in first row 

%% LDA accuracy rank - univariate
load('plmed.mat')
% load('ismale.mat')
load('ccsample.mat')
load('ccf.mat')

plmusek = plmed;
plmusek(isnan(plmusek))=0;
plmusek = plmusek(ccsample,ccf);
cfmat=zscore(plmusek);
y0=ccm;
crs=[];rrs=[];
W = (ccm)*56/43+~ccm;
for i = 1:size(cfmat,2)
    yerrs=[];randyerrs=[];ypreds=[];
    x00=cfmat(:,i);
% leave one out validation
parfor kcv=1:length(y0)
    p=[1:kcv-1,kcv+1:length(y0),kcv];
    %training
    trainnum=length(p)-1;
    trainid=p(1:trainnum);
    x=x00(trainid,:);
    y=y0(trainid);
    w0=W(trainid);
    randy = y(randperm(length(y)));
    %test
    testid=p(trainnum+1:end); 
    xt=x00(testid,:);
    yt=y0(testid);
    % Generalized Linear Model Classifier
    b = glmfit(x, y, 'binomial', 'link', 'logit','weights',w0);
    ypred = glmval(b, xt, 'logit');
    yerr = round(abs(ypred-yt));
    % Discriminant analysis classifer
    % ypred = classify(xt,x,y);
    % yerr=abs(ypred-yt);
    % Random model as bottom line 
    randypred = classify(xt,x,randy);
    randyerr=abs(randypred-yt);
    % append errors
    randyerrs=[randyerrs;randyerr];yerrs=[yerrs;yerr];
end 
% calculate correct rate
cr = 1-sum(yerrs)/length(y0);
rr = 1-sum(randyerrs)/length(y0);
crs=[crs;cr];
rrs=[rrs;rr];
end
disp(crs)
% crsR=crs;
figure;plot(Z(0.57>crs),crs(0.57>crs)*100,'ko');hold on;plot(Z(0.57<crs),crs(0.57<crs)*100,'go');bjff3;ylim([50,70]);xlim([0,4.5])

figure;plot(1:size(cfmat,2),crs*100,'ko-');bjff3;ylim([45,80]);xlim([0,115]);
%sort by accuracy and show table with feature name
[B,I] = sort(crsW,'descend');
flbl2 = flbl(:,ccf);
flbl3 = flbl2(:,I);
flbl4 = [flbl3;num2cell(B')]';


%% permutation sensitivity test 
% ageuse = pli(ccsample,cage);
% 60+
% ageuse2(ageuse<60)=0;
% ageuse2(ageuse>=60)=1;
% x0=x0(logical(ageuse2),:);
% % xg0=xg0(logical(ageuse2),:);
% y0=y0(logical(ageuse2));
maess=[];
errmats=[];
for i=1:10
    errmat=[];ypredall=[];
%     yts=[];ypreds=[];mals=[];
    parfor kcv=1:length(y0)
        p=[1:kcv-1,kcv+1:length(y0),kcv];
        trainnum=length(p)-1;
        trainid=p(1:trainnum);
        testid=p(trainnum+1:end);
        x=x0(trainid,:);
%         xg=xg0(trainid,:);
        y=y0(trainid);
%         randx = x(randperm(length(x(:,1))),:);
%         randxg = xg(randperm(length(x(:,1))),:);
        randy = y(randperm(length(y)));
        xt=x0(testid,:);
%         xgt=xg0(testid,:);
        yt=y0(testid);

%         mals=[mals;ismalignant2(testid)];
         
    %     randxt = x0(randsample(p,length(yt)),:);
        % glm
        mdl = fitglm(x,y);
        ypred = predict(mdl,xt);
        yerr=(ypred-yt);
%         mdlg = fitglm(xg,y);
%         ypredg = predict(mdlg,xgt);
%         yerrg=(ypredg-yt);
        randmdl = fitglm(x,randy);
        randypred = predict(randmdl,xt);
        randyerr=(randypred-yt);
%         randmdlg = fitglm(xg,randy);
%         randypredg = predict(randmdlg,xgt);
%         randyerrg=(randypredg-yt);
        % lasso
        lambda = 1e-03;
        [B, FitInfo] = lasso(x,y,'Alpha',0.5,'CV',10);
        idxLambda1SE = FitInfo.Index1SE;
        idxLambda1SE = FitInfo.IndexMinMSE;
        coef = B(:,idxLambda1SE);
        coef0 = FitInfo.Intercept(idxLambda1SE);
        ypred2 = xt*coef + coef0;
        yerr2= (ypred2-yt);
%         [B, FitInfo] = lasso(xg,y,'Alpha',0.5,'CV',10);
%         idxLambda1SE = FitInfo.Index1SE;
%         idxLambda1SE = FitInfo.IndexMinMSE;
%         coef = B(:,idxLambda1SE);
%         coef0 = FitInfo.Intercept(idxLambda1SE);
%         ypredg2 = xgt*coef + coef0;
%         yerrg2= (ypredg2-yt);
        [B, FitInfo] = lasso(x,randy,'Alpha',0.5,'CV',10);
        idxLambda1SE = FitInfo.Index1SE;
        idxLambda1SE = FitInfo.IndexMinMSE;
        coef = B(:,idxLambda1SE);
        coef0 = FitInfo.Intercept(idxLambda1SE);
        randypred2 = xt*coef + coef0;
        randyerr2= (randypred2-yt);
%         [B, FitInfo] = lasso(xg,randy,'Alpha',0.5,'CV',10);
%         idxLambda1SE = FitInfo.Index1SE;
%         idxLambda1SE = FitInfo.IndexMinMSE;
%         coef = B(:,idxLambda1SE);
%         coef0 = FitInfo.Intercept(idxLambda1SE);
%         randypredg2 = xgt*coef + coef0;
%         randyerrg2= (randypredg2-yt);
        % svm
        svmmdl = fitrsvm(x,y,'KernelFunction','linear','Epsilon',0.1,'Alpha',ones(size(y,1),1)*0.001 );
        ypred3 = predict(svmmdl,xt);
        yerr3=(ypred3-yt);
%         svmmdlg = fitrsvm(xg,y,'KernelFunction','linear','Epsilon',0.1,'Alpha',ones(size(y,1),1)*0.001 );
%         ypredg3 = predict(svmmdlg,xgt);
%         yerrg3=(ypredg3-yt);
        randsvmmdl = fitrsvm(x,randy,'KernelFunction','linear','Epsilon',0.1,'Alpha',ones(size(y,1),1)*0.001 );
        randypred3 = predict(randsvmmdl,xt);
        randyerr3=(randypred3-yt);
%         randsvmmdlg = fitrsvm(xg,randy,'KernelFunction','linear','Epsilon',0.1,'Alpha',ones(size(y,1),1)*0.001 );
%         randypredg3 = predict(randsvmmdlg,xgt);
%         randyerrg3=(randypredg3-yt);
        %kernel
        kernelmdl = fitrkernel(x,y,'Learner','svm','NumExpansionDimensions',2048,'KernelScale',23.8858,'Lambda',2.9317e-04,'Epsilon',2.3795);
        randMdl = fitrkernel(x,randy,'Learner','svm','NumExpansionDimensions',2048,'KernelScale',23.8858,'Lambda',2.9317e-04,'Epsilon',2.3795);
        ypred4 = predict(kernelmdl,xt);
        randypred4 = predict(randMdl,xt);
        yerr4=mean(abs(ypred4-yt));
        randyerr4=mean(abs(randypred4-yt));
    %     yerrall=[yerr(:) yerr2(:) yerr3(:)];
%         ypredall=[ypredall;[testid(:) yt(:) ypred(:) ypred2(:) ypred3(:) ypred4(:) ypredg(:) ypredg2(:) ypredg3(:) randypred(:) randypred2(:) randypred3(:) randypred4(:) randypredg(:) randypredg2(:) randypredg3(:)]];
%         errmat=[errmat;[testid(:) yerr yerr2 yerr3 yerr4 yerrg yerrg2 yerrg3 randyerr randyerr2 randyerr3 randyerr4 randyerrg randyerrg2 randyerrg3]];
        errmat=[errmat;[testid(:) yerr randyerr yerr2 randyerr2 yerr3 randyerr3 yerr4 randyerr4]];
    %     pmres=[pmres;[mean(abs(yerrall)) std(yerrall) max(abs(yerrall))]];
%         yts=[yts;yt];ypreds=[ypreds;ypred3];
    end
%     figure;plot(yts(mals),ypreds(mals),'r.');hold on;plot(yts(~mals),ypreds(~mals),'b.');plot(1:100,1:100,'k--');xticks(0:20:100);yticks(0:20:100);bjff3;xtickangle(0);
%     violin_dots({errmat(mals,6),errmat(~mals,6)},'facecolor',[.8,.8,.8;.8,.8,.8]);bjff3
    maes = mean(abs(errmat(:,2:end)));
    randcorr =[corr(errmat(:,2),errmat(:,3),'type','pearson'),corr(errmat(:,3),errmat(:,4),'type','pearson'),corr(errmat(:,6),errmat(:,7),'type','pearson'),corr(errmat(:,8),errmat(:,9),'type','pearson')];
%     disp([maes,randcorr])
    errmats=[errmats;errmat];
    maess = [maess;[maes,randcorr]];
    %glm, lasso, svm, kernel
end

maes = [[1:10]',maess(:,1:2),maess(:,9),maess(:,3:4),maess(:,10),maess(:,5:6),maess(:,11),maess(:,7:8),maess(:,12)];
maes = [maes;mean(maes)];
maes_header = {'iterations','GLM','GLM_perm','GLM_corr','Lasso','Lasso_perm','Lasso_corr','SVM','SVM_perm','SVM_corr','Kernel','Kernel_perm','Kernel_corr'};
maes2 = [maes_header;num2cell(round(maes,2))];
xlsname_permute=['permutation_sensitivity_',runtag,'.xlsx'];
writecell(maes2,xlsname_permute);

%% multivariate and gender residual comparison
load('ismale.mat')
runtag='residual check';
maess=[];errmats=[];
x00=score;
i=32;
ccm=ismale(ccsample);
x0=x00(:,1:i);
errmat=[];ypredall=[];
yts=[];ypreds=[];ccms=[];
parfor kcv=1:length(y0)
    p=[1:kcv-1,kcv+1:length(y0),kcv];
    trainnum=length(p)-1;
    trainid=p(1:trainnum);
    testid=p(trainnum+1:end);
    x=x0(trainid,:);
    y=y0(trainid);
    randy = y(randperm(length(y)));
    xt=x0(testid,:);
    yt=y0(testid);
    % glm
    mdl = fitglm(x,y);
    ypred = predict(mdl,xt);
    yerr=(ypred-yt);
    randmdl = fitglm(x,randy);
    randypred = predict(randmdl,xt);
    randyerr=(randypred-yt);
    % lasso
    lambda = 1e-03;
    [B, FitInfo] = lasso(x,y,'Alpha',0.5,'CV',10);
    idxLambda1SE = FitInfo.Index1SE;
    idxLambda1SE = FitInfo.IndexMinMSE;
    coef = B(:,idxLambda1SE);
    coef0 = FitInfo.Intercept(idxLambda1SE);
    ypred2 = xt*coef + coef0;
    yerr2= (ypred2-yt);
    [B, FitInfo] = lasso(x,randy,'Alpha',0.5,'CV',10);
    idxLambda1SE = FitInfo.Index1SE;
    idxLambda1SE = FitInfo.IndexMinMSE;
    coef = B(:,idxLambda1SE);
    coef0 = FitInfo.Intercept(idxLambda1SE);
    randypred2 = xt*coef + coef0;
    randyerr2= (randypred2-yt);
    % svm
    svmmdl = fitrsvm(x,y,'KernelFunction','linear','Epsilon',0.1,'Alpha',ones(size(y,1),1)*0.001 );
    ypred3 = predict(svmmdl,xt);
    yerr3=(ypred3-yt);
    randsvmmdl = fitrsvm(x,randy,'KernelFunction','linear','Epsilon',0.1,'Alpha',ones(size(y,1),1)*0.001 );
    randypred3 = predict(randsvmmdl,xt);
    randyerr3=(randypred3-yt);
    %kernel
    kernelmdl = fitrkernel(x,y,'Learner','svm','NumExpansionDimensions',2048,'KernelScale',23.8858,'Lambda',2.9317e-04,'Epsilon',2.3795);
    randMdl = fitrkernel(x,randy,'Learner','svm','NumExpansionDimensions',2048,'KernelScale',23.8858,'Lambda',2.9317e-04,'Epsilon',2.3795);
    ypred4 = predict(kernelmdl,xt);
    randypred4 = predict(randMdl,xt);
    yerr4=mean(abs(ypred4-yt));
    randyerr4=mean(abs(randypred4-yt));
    errmat=[errmat;[testid(:) yerr randyerr yerr2 randyerr2 yerr3 randyerr3 yerr4 randyerr4]];
    yts=[yts;yt];ypreds=[ypreds;ypred3];ccms=[ccms;ccm(testid)];
end
maes = mean(abs(errmat(:,2:end)));
randcorr =[corr(errmat(:,2),errmat(:,3),'type','pearson'),corr(errmat(:,3),errmat(:,4),'type','pearson'),corr(errmat(:,6),errmat(:,7),'type','pearson'),corr(errmat(:,8),errmat(:,9),'type','pearson')];
errmats=[errmats;errmat];
maess = [maess;[maes,randcorr]];

maes = [1,maess(:,1:2),maess(:,9),maess(:,3:4),maess(:,10),maess(:,5:6),maess(:,11),maess(:,7:8),maess(:,12)];
% maes = [maes;mean(maes)];
maes_header = {'iterations','GLM','GLM_perm','GLM_corr','Lasso','Lasso_perm','Lasso_corr','SVM','SVM_perm','SVM_corr','Kernel','Kernel_perm','Kernel_corr'};
maes2 = [maes_header;num2cell(round(maes,2))];
% xlsname_permute=['Npc_vs_MAE_',runtag,'.xlsx'];
% writecell(maes2,xlsname_permute);
ccms=logical(ccms);
plot(yts,ypreds,'k.','markersize',10);hold on;plot(1:100,1:100,'k--');bjff3;
figure;plot(yts(ccms),ypreds(ccms),'b.','markersize',10);hold on;plot(yts(~ccms),ypreds(~ccms),'r.','markersize',10);plot(1:100,1:100,'k--');bjff3;
residual = errmat(:,6);
figure;plot(yts(ccms),residual(ccms),'b.','markersize',10);hold on;plot(yts(~ccms),residual(~ccms),'r.','markersize',10);plot(1:100,zeros(1,100),'k--');bjff3;

%% build gender-independent feature model using male only and female only population to test sample size effect
load('isback.mat')
load('ismale.mat')
load('iswhite.mat')
load('ccfeat_global.mat')
load('ccfeat_female.mat')
load('ccfeat_male.mat')
load('age.mat')
load('plmed.mat')

ccsample=isback&iswhite;
maess=[];errmats=[];%33 for white back; 34 for multisec

%% male population, global feature, remove nan sample, manual pc
ccfeat=cf;
x0=plmed(ccsample&ismale,ccfeat);y0=age(ccsample&ismale);
dopca=1; removenansample=1; replacenanvalue=0;
if removenansample;[row, ~] = find(isnan(x0));delidx=unique(row);x0(delidx,:)=[];y0(delidx)=[];end
if replacenanvalue;x0(isnan(x0))=0;end
if dopca;x0=zscore(x0);[coeff,score,latent]=pca(x0);vcum=cumsum(latent)/sum(latent);v95num=find(vcum-0.95>0,1,'first');end 
% x_reduced=score(:,1:v95num);
disp(size(x0,2))
disp(length(y0))
disp(sum(v95num)) %display 
pc_search = [v95num-2,v95num-1,v95num,v95num+1,v95num+2];
maesspc=[];errmatspc=[];
for i = 1:length(pc_search)
    x_reduced=score(:,1:pc_search(i));
    [es,randcorr,errmat,predictglm,predictlas,predictsvm,predictker] = agingmodel(x_reduced,y0);
%     figure;plot(predictsvm(:,2),predictsvm(:,1),'k.');hold on;plot(1:100,1:100,'k--');bjff3
    errmatspc=[errmatspc;errmat];
    maesspc = [maesspc;[maes,randcorr]];
end
mae = [[1:size(maesspc,1)]',maesspc(:,1:2),maesspc(:,9),maesspc(:,3:4),maesspc(:,10),maesspc(:,5:6),maesspc(:,11),maesspc(:,7:8),maesspc(:,12)];

%% female, global, removenansample, manual pc
ccfeat=cf;
x0=plmed(ccsample&ismale,ccfeat);y0=age(ccsample&ismale);
dopca=1; removenansample=1; replacenanvalue=0;
if removenansample;[row, ~] = find(isnan(x0));delidx=unique(row);x0(delidx,:)=[];y0(delidx)=[];end
if replacenanvalue;x0(isnan(x0))=0;end
if dopca;x0=zscore(x0);[coeff,score,latent]=pca(x0);vcum=cumsum(latent)/sum(latent);v95num=find(vcum-0.95>0,1,'first');end 
% x_reduced=score(:,1:v95num);
disp(length(y0))
disp(size(x0,2))
disp(sum(v95num)) %display 
pc_search = [v95num-2,v95num-1,v95num,v95num+1,v95num+2];
maesspc=[];errmatspc=[];
for i = 1:length(pc_search)
    x_reduced=score(:,1:pc_search(i));
    [maes,randcorr,errmat,predictglm,predictlas,predictsvm,predictker] = agingmodel(x_reduced,y0);
%     figure;plot(predictsvm(:,2),predictsvm(:,1),'k.');hold on;plot(1:100,1:100,'k--');bjff3
    errmatspc=[errmatspc;errmat];
    maesspc = [maesspc;[maes,randcorr]];
end
mae = [[1:size(maesspc,1)]',maesspc(:,1:2),maesspc(:,9),maesspc(:,3:4),maesspc(:,10),maesspc(:,5:6),maesspc(:,11),maesspc(:,7:8),maesspc(:,12)];

%% all
ccfeat=cf;
x0=plmed(ccsample&age>60,ccfeat);y0=age(ccsample&age>60);
dopca=1; removenansample=0; replacenanvalue=1;
if removenansample;[row, ~] = find(isnan(x0));delidx=unique(row);x0(delidx,:)=[];y0(delidx)=[];end
if replacenanvalue;x0(isnan(x0))=0;end
if dopca;x0=zscore(x0);[coeff,score,latent]=pca(x0);vcum=cumsum(latent)/sum(latent);v95num=find(vcum-0.95>0,1,'first');end 
% x_reduced=score(:,1:v95num); for i = 1:length(pc_search)
    x_reduced=score(:,1:pc_search(i));
    [maes,randcorr,errmat,predictglm,predictlas,predictsvm,predictker] = agingmodel(x_reduced,y0);
%     figure;plot(predictsvm(:,2),predictsvm(:,1),'k.');hold on;plot(1:100,1:100,'k--');bjff3
    errmatspc=[errmatspc;errmat];
    maesspc = [maesspc;[maes,randcorr]];
% end
mae = [[1:size(maesspc,1)]',maesspc(:,1:2),maesspc(:,9),maesspc(:,3:4),maesspc(:,10),maesspc(:,5:6),maesspc(:,11),maesspc(:,7:8),maesspc(:,12)];


% maes = [maes;mean(maes)];
mae_header = {'iterations','GLM','GLM_perm','GLM_corr','Lasso','Lasso_perm','Lasso_corr','SVM','SVM_perm','SVM_corr','Kernel','Kernel_perm','Kernel_corr'};
mae2 = [mae_header;num2cell(round(mae,2))];
xlsname_permute=['gender',runtag,'.xlsx'];
writecell(mae2,xlsname_permute);

plot(maes(:,1),maes(:,2),'-');hold on; plot(maes(:,1),maes(:,5),'-'); plot(maes(:,1),maes(:,8),'-');plot(maes(:,1),maes(:,11),'-');bjff3
legend('GLM','Lasso','SVM','Kernel')

%% rank feature to differentiate male vs female
ccsample = iswhite&isback;
BC=ismale(ccsample);
X=plmed(ccsample,:)';
cps=zeros(1,1090);
Nsample=50;
I = rankfeatures(X,BC,NumberOfIndices=Nsample);
for i = 1:Nsample
    try
    C = classify(X(I(1:i),:)',X(I(1:i),:)',double(BC));
    catch
        disp(i)
        cps(i) = i;
    end
    cp = classperf(BC,C);
    cps(i) = cp.CorrectRate;
end
plot(1:Nsample,cps(1:Nsample),'k+')

%%
ps=ones(1,108);
for i=1:108
[h,p]=ttest2(cfmat(ccm,i),cfmat(~ccm,i));
ps(i)=p;
end

