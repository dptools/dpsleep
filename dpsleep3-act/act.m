function act(read_dir, output_dir,outdp_dir study, subject, ref_date)

%% Study and subject 
stdy=study; 
sb1=subject; 
mn=datenum(ref_date);   % Choose 1st day of the study as the actigraphy reference

%% Input and output directory
display('New Pipeline.')
display('Updated Pipeline.');
display('Checking input directory.');
adr=read_dir;
% Check if the path is properly formatted
if ~ endsWith(adr, '/')
    adr = strcat(adr, '/');
end

display('Checking output directory.');
adrq=output_dir;
% Check if the path is properly formatted
if ~ endsWith(adrq, '/')
    adrq = strcat(adrq, '/');
end

%% Find the files related to the day + some days before and after that day
display('Finding files.');
d3=dir(strcat(adr,'*.mat'));

files_len = length(d3);
% Exit if there are no files to read
if files_len == 0
    display('Files do not exist under this directory.');
    exit(1);
end

%% Parameters
display('Initializing parameters.');
w1=100;              % Movig window size (mins) for sleep
w2=15;               % Movig window size (mins) for fractions
w=10;                % Adjust beginning and end of sleep
%wdb=5;   wdf=223;   % Window of Analysis
fth=0.3;             % Frequency threshhold to calculate mean for sleep  
wsd1=150;            % Moving window size to find watch off long
wsd2=60;             % Moving window size to find watch off short
thra=0.05;          % Threshhold to find active/non-active minutes  
fcth=1.2;             % Frequency threshhold to find acive minutes 
thrw=.05; thrr=.05;    %Thresholds for walking and running

d4=extractfield(d3,'name');
d5=split(d4,'.'); d61=d5(:,:,1);
d51=split(d61,'_'); d6=d51(:,:,2);
d7=datenum(d6);
[d7s,ids]=sort(d7,'ascend');
d7sr=d7s-mn+1;
wdf=max(d7sr); wdb=0;
d4s=d4(ids);

indx=3*ones(1440*(wdb+wdf),1);  indb=zeros(1440*(wdb+wdf),1); indxp=zeros(1440*(wdb+wdf),1); indxbp=zeros(1440*(wdb+wdf),1); indxpd=zeros(1440*(wdb+wdf),1); 
indxx=3*ones(1440*(wdb+wdf),1); inds=zeros(1440*(wdb+wdf),1);  indx_slp=2*ones(1440*(wdb+wdf),1); indx_phil=2*ones(1440*(wdb+wdf),1); ind90_slp=2*ones(1440*(wdb+wdf),1);
  indpf=ones(1440*(wdb+wdf),1);

%% Extract data both time and frequency
td=(1/1440):(1/1440):(wdb+wdf);
tdd=(1/1440):(1/1440):(wdb+wdf);
idar=d7sr;
ida=1:wdf;
[idall,ia,ib]=intersect(ida,idar);
mxtal1=[]; mytal1=[]; mztal1=[]; sxtal1=[]; sytal1=[]; sztal1=[]; bff=[]; mgf=[]; spf=[];
i2=0;
for i1=1:length(ida)
    if ~ismember(i1,ia)       
       bff0=zeros(1440,1);  spf0=zeros(1440,1);
       mxt=zeros(1,1440); myt=zeros(1,1440); mzt=zeros(1,1440); mgt=zeros(1,1440);
       sxt=zeros(1,1440); syt=zeros(1,1440); szt=zeros(1,1440);
    else
        i2=i2+1;
        i=ib(i2);
        load(strcat(adr,d4s{1,i}))                    
        %% Mean spectrum       
        mfn=f_all(1,:);
        px_all(:,mfn<fth)=NaN;  py_all(:,mfn<fth)=NaN;  pz_all(:,mfn<fth)=NaN;
        xmf=nanmean(px_all,2); % Read spectrum
        ymf=nanmean(py_all,2);
        zmf=nanmean(pz_all,2);
        spf0=(xmf.^2+ymf.^2+zmf.^2).^(.5);
        bff0=btm_all; %#ok<*AGROW>                
        mxf=round(max(f_all,[],2)*4)/2;   % Calculate the sampling frequency at each minute
        %% Standard deviation
        xt=reshape(axyz_all(:,1),mxf(1)*60,1440);
        yt=reshape(axyz_all(:,2),mxf(1)*60,1440);
        zt=reshape(axyz_all(:,3),mxf(1)*60,1440);
        lgt=reshape(light_all,mxf(1)*60,1440);
        mxt=nanmean(xt); myt=nanmean(yt);  mzt=nanmean(zt);  mgt=nanmean(lgt);
        sxt=nanstd(xt); syt=nanstd(yt);  szt=nanstd(zt);        
    end
    mxtal1=[mxtal1;mxt']; mytal1=[mytal1;myt']; mztal1=[mztal1;mzt'];  mgf=[mgf;mgt'];
    sxtal1=[sxtal1;sxt']; sytal1=[sytal1;syt']; sztal1=[sztal1;szt'];
    spf=[spf;spf0]; bff=[bff;bff0];
end
 bf=bff; 

%% Find watch off epochs
sxtal=sxtal1;   sytal=sytal1; sztal=sztal1;
sxtal(isnan(sxtal))=0; sytal(isnan(sytal))=0; sztal(isnan(sztal))=0;
sxmv1=movmean(sxtal,[0 wsd1]);
symv1=movmean(sytal,[0 wsd1]);
szmv1=movmean(sztal,[0 wsd1]);
sxmb1=movmean(sxtal,[wsd1 0]);
symb1=movmean(sytal,[wsd1 0]);
szmb1=movmean(sztal,[wsd1 0]);
smv1=(sxmv1.^2+symv1.^2+szmv1.^2).^(.5);
smb1=(sxmb1.^2+symb1.^2+szmb1.^2).^(.5);
%% Find moving average of spectrum (forward) Short window
saf=movmean(spf,[0 w2]);
smt=0.018525;
saf1=saf;
saf1(smv1<smt | smb1<smt)=0;      % watch off long
idz=find(saf1>0.0000);       % Delete zero values
safz=saf1(idz);   tafz=td(idz);    

%% Find moving average of spectrum (backward) Short window
sab=movmean(spf,[w2 0]);  

%% Find low activity modes
prbf2=prctile(safz,10);
prbf3=prctile(safz,25);
prbf5=prbf3;
prbf7=prctile(safz,50);
prbf8=prctile(safz,75);

%% Find moving average of spectrum 1min window
saff=movmean(spf,[0 1]);
saff1=saff;
saff1(smv1<smt | smb1<smt)=0;      % watch off long
idz=find(saff1>0.000);       % Delete zero values
saffz=saff1(idz);   taffz=td(idz);

%% Find low activity modes
prbff1=prctile(saffz,5);
prbff2=prctile(saffz,10);
prbff3=prctile(saffz,15);
prbff4=prctile(saffz,20);
prbff5=prctile(saffz,25);
prbff6=prctile(saffz,30);
prbff7=prctile(saffz,35);
prbff8=prctile(saffz,40);
prbff9=prctile(saffz,45);
prbff10=prctile(saffz,50);
prbff11=prctile(saffz,55);
prbff12=prctile(saffz,60);
prbff13=prctile(saffz,65);
prbff14=prctile(saffz,70);
prbff15=prctile(saffz,75);
prbff16=prctile(saffz,80);
prbff17=prctile(saffz,85);
prbff18=prctile(saffz,90);
prbff19=prctile(saffz,95);

%% Watch off long
indx(smv1<smt | smb1<smt)=0;      % watch off long
indxx(smv1<smt | smb1<smt)=0;      % watch off long
%% Sleep Epoch Detection
inds1=inds; inds2=inds;

%% Activity Scores 
inds((saf<prbf5 | sab<prbf5))=1;  % <25%
inds((saf<prbf2 | sab<prbf2))=1;  % <10%
inds((saf>prbf7 | sab>prbf7))=-.75; % >50%
inds((saf>prbf8 | sab>prbf8))=-1;   % >75%
inds(smv1<smt | smb1<smt)=0;        % watch off long

%% Low Activity Candidates 100mins and 60mins windows
isf=movmean(inds,[0 100]);   isb=movmean(inds,[100 0]);
isfs=movmean(inds,[0 60]);   isbs=movmean(inds,[60 0]);
zisf=isf(isf>0);   
inds1((isf>=0 & isbs>=0 ) |( isb>=0 & isfs>=0) )=1;
ind90_slp(inds1==1)=1;
ind90_slp(smv1<smt | smb1<smt)=0;  % watch off long
inds1(smv1<smt | smb1<smt)=0;      % watch off long
ind90_slp(inds1==1)=1;

%% Sleep Candidates 90mins window
issf=movmean(inds1,[0 90]);   issb=movmean(inds1,[90 0]); 
inds1((issf<=.75) |(issb<=.75) )=0;
inds1((issf>.75) |(issb>.75) )=1;
indx_slp(inds1==1)=1;
indx_slp(smv1<smt | smb1<smt)=0;      % watch off long

%% Shift to 6pm-6pm
bgn=18*60;  endt=6*60;
indx_slp2=[zeros(1440,1);indx_slp];
indx_slp1=indx_slp2(bgn+1:end-endt);
indf_slp=reshape(indx_slp1,1440,wdb+wdf);
ind90_slp2=[zeros(1440,1);ind90_slp];
ind90_slp1=ind90_slp2(bgn+1:end-endt);
indf90_slp=reshape(ind90_slp1,1440,wdb+wdf);
indx(inds1==1)=6;
%% Reshape Light
mgf(isnan(mgf))=0;
mgf2=[zeros(1440,1);mgf];
mgf3=mgf2(bgn+1:end-endt);
indf_mgf=reshape(mgf3,1440,wdb+wdf);
%% Color codes
%% 50% Highest Activity 
indx((saf>prbf7 | sab>prbf7) )=4;
indxx((saf>prbf7 | sab>prbf7) )=4;
%% 75% Highest Activity 
indx((saf>prbf8 | sab>prbf8) )=5;
indxx((saf>prbf8 | sab>prbf8) )=5;
%% 25% Lowest Activity 
indx((saf<=prbf3 | sab<=prbf3) )=2;
indx(smv1<smt | smb1<smt)=0;      % watch off long
indxx((saf<=prbf3 | sab<=prbf3) )=2;
indxx(smv1<smt | smb1<smt)=0;      % watch off long
%% 10% Lowest Activity 
indx((saf<=prbf2 | sab<=prbf2))=1;
indx(smv1<smt | smb1<smt)=0;      % watch off long
indxx((saf<=prbf2 | sab<=prbf2))=1;
indxx(smv1<smt | smb1<smt)=0;      % watch off long

%% Mark button press
bfw=movmean(bf,[0 3]);
indx(bfw>0 & bfw<500)=7;
indx(smv1<smt | smb1<smt)=0;
indxx(bfw>0 & bfw<500)=6;
indx(smv1<smt | smb1<smt)=0;
%% Find relative dates
t0=mn-datenum(2015,1,31);
drl=floor(mn)-mn+1;
drls=drl-wdb;   drlp=drl+wdf;
indx2=[zeros(1440,1);indx];
indx1=indx2(bgn+1:end-endt);
indxx2=[zeros(1440,1);indxx];
indxx1=indxx2(bgn+1:end-endt);

%%  Activity Levels
indxp(saff>prbff19)=20;
indxp(saff<=prbff19)=19;
indxp(saff<=prbff18)=18;
indxp(saff<=prbff17)=17;
indxp(saff<=prbff16)=16;
indxp(saff<=prbff15)=15;
indxp(saff<=prbff14)=14;
indxp(saff<=prbff13)=13;
indxp(saff<=prbff12)=12;
indxp(saff<=prbff11)=11;
indxp(saff<=prbff10)=10;
indxp(saff<=prbff9)=9;
indxp(saff<=prbff8)=8;
indxp(saff<=prbff7)=7;
indxp(saff<=prbff6)=6;
indxp(saff<=prbff5)=5;
indxp(saff<=prbff4)=4;
indxp(saff<=prbff3)=3;
indxp(saff<=prbff2)=2;
indxp(saff<=prbff1)=1;
%indxp(smv1<smt | smb1<smt)=0;      % watch off long
%%
safff=saff;
safff(safff>=prbff15)=prbff15;
indxpd=100*safff/(max(safff(:)));

%% Mark button press
bfw1=movmean(bf,[0 1]);
indxbp(bfw1>0)=1;
%% Daily format
indxp2=[zeros(1440,1);indxp];
indxp1=indxp2(bgn+1:end-endt);
indf_xp=reshape(indxp1,1440,wdb+wdf);
indxbp2=[zeros(1440,1);indxbp];
indxbp1=indxbp2(bgn+1:end-endt);
indf_bp=reshape(indxbp1,1440,wdb+wdf);
indxpd2=[zeros(1440,1);indxpd];
indxpd1=indxpd2(bgn+1:end-endt);
indf_xpd=reshape(indxpd1,1440,wdb+wdf);
indsaf2=[zeros(1440,1);saff];
indsaf1=indsaf2(bgn+1:end-endt);
indf_saf=reshape(indsaf1,1440,wdb+wdf);
indf_act=reshape(indxx1,1440,wdb+wdf);
indp_act=reshape(indxx,1440,wdb+wdf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DPdash Scores %%%%%%%%%%%%%%%%%%%%%%
indf_score=indp_act;
ndy=length(indf_score(1,:));
indv_scr=reshape(indf_score,1,ndy*1440);
indv_scr(indv_scr==0)=NaN; 
indf_scr=reshape(indv_scr,1440,ndy);
indd_score=nanmean(indf_scr)';
indf_scr2=reshape(indv_scr,60,24*ndy);
indh_score1=nanmean(indf_scr2);
if1=find(~isnan(indh_score1),1,'first');
ife=24*ceil(if1/24);
ifs=ife-if1+1;
indh_ppy1=[NaN(1,ifs) indh_score1(1:end-ifs)];
indh_score=reshape(indh_score1,24,ndy)';
indh_ppy=reshape(indh_ppy1,24,ndy)';
indh_pp=indh_ppy';
indh_pp(~isnan(indh_pp))=1;
indd_pp=nansum(indh_pp);
indd_ppy1=zeros(size(indd_pp));
indd_ppy1(indd_pp>23)=1;
il1=find(~isnan(indd_pp),1,'last');
if indd_pp(il1)>12
    indd_ppy1(il1)=1;
end
indd_ppy=[indd_ppy1;cumsum(indd_ppy1,'omitnan')];
indd_ppy=indd_ppy';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DPdash format. %%%%%%%%%%%%%%%%%%%%%%%%
dys=[1:1:ndy]';
mdys=dys+mn-1;
wdys=weekday(mdys);
jdys=(mdys-datenum('1970-01-01'))*1000*3600*24;
day=num2str(dys);
reftime=num2str(jdys);
weekday1=num2str(wdys);
timeofday=[];
for k=1:ndy
    timeofday=[timeofday;'00:00:00'];
end
tab_ddp=table(reftime,day,timeofday,weekday1);
try
   % tab_ddp.Properties.VariableNames{'days'} = 'day';
    tab_ddp.Properties.VariableNames{'weekday1'} = 'weekday';
catch ME

end
score=indh_score;
score2=round(100*score)/100;
score1=string(score2);
tab_scr_hr1 =array2table(score1);
tab_scr_hr1.Properties.VariableNames={'activityScore_hour01','activityScore_hour02', 'activityScore_hour03','activityScore_hour04','activityScore_hour05','activityScore_hour06','activityScore_hour07','activityScore_hour08','activityScore_hour09','activityScore_hour10','activityScore_hour11','activityScore_hour12','activityScore_hour13','activityScore_hour14','activityScore_hour15','activityScore_hour16','activityScore_hour17','activityScore_hour18','activityScore_hour19','activityScore_hour20','activityScore_hour21','activityScore_hour22','activityScore_hour23','activityScore_hour24'};
tab_scr_hr=[tab_ddp tab_scr_hr1];
%%
display('Checking output directory for DPDash.');
adrq0=outdp_dir;
% Check if the path is properly formatted
if ~ endsWith(adrq0, '/')
    adrq0 = strcat(adrq0, '/');
end
adrq00=strcat(adrq0,'accel');
if exist(adrq00,'dir')~=7
    mkdir(adrq00)
end
writetable(tab_scr_hr,strcat(adrq00,'/',study,'-',sb1,'-actigraphy_GENEActiv_accel_activityScores_hourly-day1to',num2str(ndy),'.csv'),'Delimiter',',','QuoteStrings',false)
%% 
scoref=indh_ppy;
score4=round(100*scoref)/100;
score3=string(score4);
tab_scr_hr2 =array2table(score3);
tab_scr_hr2.Properties.VariableNames={'activityScore_hour01','activityScore_hour02', 'activityScore_hour03','activityScore_hour04','activityScore_hour05','activityScore_hour06','activityScore_hour07','activityScore_hour08','activityScore_hour09','activityScore_hour10','activityScore_hour11','activityScore_hour12','activityScore_hour13','activityScore_hour14','activityScore_hour15','activityScore_hour16','activityScore_hour17','activityScore_hour18','activityScore_hour19','activityScore_hour20','activityScore_hour21','activityScore_hour22','activityScore_hour23','activityScore_hour24'};
tab_scr_hrf=[tab_ddp tab_scr_hr2];
writetable(tab_scr_hrf,strcat(adrq00,'/',study,'-',sb1,'-actigraphy_GENEActiv_accshft_activityScores_hourly-day1to',num2str(ndy),'.csv'),'Delimiter',',','QuoteStrings',false)
%% 
score5=indd_ppy;
scorep=string(score5);
tab_scr_hr3 =array2table(scorep);
tab_scr_hr3.Properties.VariableNames={'pay','payc'};
tab_scr_hrp=[tab_ddp tab_scr_hr3];
writetable(tab_scr_hrp,strcat(adrq00,'/',study,'-',sb1,'-actigraphy_GENEActiv_pay_daily-day1to',num2str(ndy),'.csv'),'Delimiter',',','QuoteStrings',false)
%% Build data and figures folder 
outp=adrq;
if exist(outp,'dir')~=7
   mkdir(outp) 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save mat file %%%%%%%%%%%%%%%%%%%%%%%%%
save(strcat(outp,stdy,'-',sb1,'-geneactiv_mtl3.mat'),'indf_slp','indf90_slp','indf_act','indf_mgf','indf_xp','indf_bp','indf_saf','indf_xpd')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% You can comment from here to bypass the Plots (except for the last 4 rows) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot with sleep %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ff1=figure(1);
set(gcf,'position',get(0,'screensize'))
set(gcf,'color','white')
set(gcf,'defaultAxesFontName','Arial')
if isempty(find(indx1==0)) %#ok<*EFIND>
    clrs=[0.5 1 1; 1 0.5 1;.9 .9 .9; 0 0.5 0.4; 0 0 1; 1 1 .75; 0 0 0];
else 
    clrs=[1 1 1; 0.5 1 1; 1 0.5 1;.9 .9 .9; 0 0.5 0.4; 0 0 1; 1 1 .75; 0 0 0];
end
colormap(gca,clrs)
clr_in1=clrs;
indf=reshape(indx1,1440,wdb+wdf);
hhind=imagesc(indf);
axi = ancestor(hhind, 'axes');
indt=indf_slp;
indt(indt==2)=1;
set(hhind,'AlphaData',indt)
xrule = axi.XAxis;
h25fpp=colorbar;
yrule = axi.YAxis;
set(gca, 'FontWeight','b')
epc_lbl{1,1}={' 6PM',' 7PM',' 8PM',' 9PM','10PM','11PM','12AM',' 1AM',' 2AM',' 3AM',' 4AM',' 5AM',' 6AM',' 7AM',' 8AM',' 9AM','10AM','11AM','12PM',' 1PM',' 2PM',' 3PM',' 4PM',' 5PM',' 6PM'};
epochm={[0 24];[9 17];[18 22];[22 2];[2 6]};
ipc=1;
ylim([epochm{ipc,1}(1,1)*60+.5 epochm{ipc,1}(1,2)*60+.5])
set(gca,'YTick',(epochm{ipc,1}(1,1)*60+0.5):60:(epochm{ipc,1}(1,2)*60+0.5),'YTickLabel',epc_lbl{ipc,1}(1:1:end))           
set(gca,'XTick',0.5:2:length(indf(1,:))+.5,'XTickLabel',drls-1:2:drlp,'XTickLabelRotation',90)
xrule.FontSize = 18;
yrule.FontSize = 18;
xlabel('Relative Days','FontWeight','b','FontSize',18)
grid on
if isempty(find(indx1==0))
    lblf2={'ActL<10%';'10%<ActL<25%';'25%<ActL<50%'; '50%<ActL<75%'; 'ActL>75%';'Sleep';'Button'};
else 
    lblf2={'Watch Off';'ActL<10%';'10%<ActL<25%';'25%<ActL<50%'; '50%<ActL<75%'; 'ActL>75%';'Sleep';'Button'};
end
set(h25fpp,'YTick',.5:.87:8.5,'YTickLabel',lblf2,'FontSize',18)
%% Save figures
outp=adrq;
if exist(outp,'dir')~=7
   mkdir(outp) 
end
savefig(strcat(outp,stdy,'-',sb1(1:3),'-geneactiv_mtl3s.fig'))
img = getframe(gcf);
imwrite(img.cdata, strcat(outp,stdy,'-',sb1(1:3),'-geneactiv_mtl3s.png'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot without sleep %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ff2=figure(2);
set(gcf,'position',get(0,'screensize'))
set(gcf,'color','white')
set(gcf,'defaultAxesFontName','Arial')
if isempty(find(indxx1==0)) %#ok<*EFIND>
    clrs=[0 0 1;0 1 1; 0 1 0; 1 0.65 0;1 0 0; 0 0 0];
else 
    clrs=[1 1 1; 0 0 1;0 1 1; 0.466 .67 0.188; 1 0.65 0;1 0 0; 1 1 0];
end
colormap(gca,clrs)
clr_in1=clrs;
indf=indf_act;
indff=flip(indf,1);
hhind=imagesc(indff);
axi = ancestor(hhind, 'axes');
xrule = axi.XAxis;
%h35fpp=colorbar;
yrule = axi.YAxis;
set(gca, 'FontWeight','b')
epc_lbl{1,1}={' 6PM',' 3PM', '12PM', ' 9AM',' 6AM',' 3AM','12AM',' 9PM',' 6PM'};
epochm={[0 24];[9 17];[18 22];[22 2];[2 6]};
ipc=1;
ylim([0*60+.5 24*60+.5])
set(gca,'YTick',(0*60+0.5):3*60:(24*60+0.5),'YTickLabel',epc_lbl{ipc,1}(1:1:end),'TickDir','out', 'Linewidth',4,'FontName','Arial','FontWeight','b','FontSize',32,'Ygrid','off')
set(gca,'XTick',[],'XTickLabel',[],'XTickLabelRotation',0,'FontSize',32)
%xlabel('Relative Days','FontWeight','b','FontSize',32,'FontName','Arial')

a = gca;                                % get handle to current axes
set(a,'box','off','color','none')           % set box property to off and remove background color
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'linewidth',4);   % create new, empty axes with box but without ticks
axes(a)          % set original axes as active
linkaxes([a b])  % link axes in case of zooming

%grid on
if isempty(find(indxx1==0))
    lblf2={'0%<ActL<10%';'10%<ActL<25%';'25%<ActL<50%'; '50%<ActL<75%'; '75%<ActL<100%';'Button'};
else 
    lblf2={'Watch Off';'0%<ActL<10%';'10%<ActL<25%';'25%<ActL<50%'; '50%<ActL<75%'; '75%<ActL<100%';'Button'};
end
%set(h35fpp,'YTick',.5:.85:7.5,'YTickLabel',lblf2,'FontSize',32)
%% Save figures
savefig(strcat(outp,stdy,'-',sb1(1:3),'-geneactiv_mtl3p.fig'))
img = getframe(gcf);
imwrite(img.cdata, strcat(outp,stdy,'-',sb1(1:3),'-geneactiv_mtl3p.png'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot without sleep %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ff4=figure(4);
set(gcf,'position',get(0,'screensize'))
set(gcf,'color','white')
set(gcf,'defaultAxesFontName','Arial')
clps=unique(indf_xp);
excl=clps;
clrex=jet(length(excl));
clrx=[1 1 1;clrex(1:end,:)];
colormap(gca,clrx)
hhind4=imagesc(indf_xp);
axi = ancestor(hhind, 'axes');
xrule = axi.XAxis;
h35fppx=colorbar;
yrule = axi.YAxis;
set(gca, 'FontWeight','b')
epc_lbl{1,1}={' 6PM',' 9PM', '12AM', ' 3AM',' 6AM',' 9AM','12PM',' 3PM',' 6PM'};
epochm={[0 24];[9 17];[18 22];[22 2];[2 6]};
ylim([0*60+.5 24*60+.5])
set(gca,'YTick',(0*60+0.5):3*60:(24*60+0.5),'YTickLabel',epc_lbl{1,1}(1:1:end),'TickDir','out', 'Linewidth',4,'FontName','Arial','FontWeight','b','FontSize',32,'Ygrid','off' )
set(gca,'XTick',[.5 2.5:2:length(indf(1,:))+.5],'XTickLabel',[1 drls+2:2:drlp],'XTickLabelRotation',0,'FontSize',32)

a = gca;                                % get handle to current axes
set(a,'box','off','color','none')           % set box property to off and remove background color
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'linewidth',4);   % create new, empty axes with box but without ticks
axes(a)          % set original axes as active
linkaxes([a b])  % link axes in case of zooming

xlabel('Relative Days','FontWeight','b','FontSize',32,'FontName','Arial')
%grid on
if isempty(find(indxx1==0))
    lblf2={'0%<ActL<10%';'10%<ActL<25%';'25%<ActL<50%'; '50%<ActL<75%'; '75%<ActL<100%';'Button'};
else 
    lblf2={'Watch Off';'0%<ActL<10%';'10%<ActL<25%';'25%<ActL<50%'; '50%<ActL<75%'; '75%<ActL<100%';'Button'};
end
%set(h35fpp,'YTick',.5:.85:7.5,'YTickLabel',lblf2,'FontSize',32)
%% Save figures
savefig(strcat(outp,stdy,'-',sb1(1:3),'-geneactiv_mtl3act.fig'))
img = getframe(gcf);
imwrite(img.cdata, strcat(outp,stdy,'-',sb1(1:3),'-geneactiv_mtl3act.png'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot sleep before 90mins window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ff33=figure(33);
set(gcf,'position',get(0,'screensize'))
set(gcf,'color','white')
clrs=[1 1 1; 0 0 1; 0.9 0.9 0.9];
colormap(gca,clrs)
clr_in1=clrs;
indf_slpp=flip(indf90_slp,1);
hhind=imagesc(indf_slpp);
axi = ancestor(hhind, 'axes');
xrule = axi.XAxis;
yrule = axi.YAxis;
set(gca, 'FontWeight','b')
epc_lbl{1,1}={' 6PM',' 3PM', '12PM', ' 9AM',' 6AM',' 3AM','12AM',' 9PM',' 6PM'};
epochm={[0 24];[9 17];[18 22];[22 2];[2 6]};
ipc=1;
ylim([epochm{ipc,1}(1,1)*60+.5 epochm{ipc,1}(1,2)*60+.5])
set(gca,'YTick',(0*60+0.5):3*60:(24*60+0.5),'YTickLabel',epc_lbl{ipc,1}(1:1:end),'TickDir','out', 'Linewidth',4,'FontName','Arial','FontWeight','b','FontSize',32,'Ygrid','off')           
set(gca,'XTick',[],'XTickLabel',[],'XTickLabelRotation',0,'FontSize',32)

a = gca;                                % get handle to current axes
set(a,'box','off','color','none')           % set box property to off and remove background color
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'linewidth',4);   % create new, empty axes with box but without ticks
axes(a)          % set original axes as active
linkaxes([a b])  % link axes in case of zooming

%xlabel('Relative Days','FontWeight','b','FontSize',32)
%grid on
lblf2={'Watch Off';'Sleep';'Awake'};
%set(h35fpp,'YTick',.3:.7:2.3,'YTickLabel',lblf2,'FontSize',32)
%% Save figures
savefig(strcat(outp,stdy,'-',sb1(1:3),'-geneactiv_mtl33ss.fig'))
img = getframe(gcf);
imwrite(img.cdata, strcat(outp,stdy,'-',sb1(1:3),'-geneactiv_mtl33ss.png'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot sleep %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ff3=figure(3);
set(gcf,'position',get(0,'screensize'))
set(gcf,'color','white')
clrs=[1 1 1; 0 0 1; 0.9 0.9 0.9];
colormap(gca,clrs)
clr_in1=clrs;
indf_slpp=flip(indf_slp,1);
hhind=imagesc(indf_slpp);
axi = ancestor(hhind, 'axes');
xrule = axi.XAxis;
yrule = axi.YAxis;
set(gca, 'FontWeight','b')
epc_lbl{1,1}={' 6PM',' 3PM', '12PM', ' 9AM',' 6AM',' 3AM','12AM',' 9PM',' 6PM'};
epochm={[0 24];[9 17];[18 22];[22 2];[2 6]};
ipc=1;
ylim([epochm{ipc,1}(1,1)*60+.5 epochm{ipc,1}(1,2)*60+.5])
set(gca,'YTick',(0*60+0.5):3*60:(24*60+0.5),'YTickLabel',epc_lbl{ipc,1}(1:1:end),'TickDir','out', 'Linewidth',4,'FontName','Arial','FontWeight','b','FontSize',32,'Ygrid','off')           
set(gca,'XTick',[],'XTickLabel',[],'XTickLabelRotation',0,'FontSize',32)

a = gca;                                % get handle to current axes
set(a,'box','off','color','none')           % set box property to off and remove background color
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'linewidth',4);   % create new, empty axes with box but without ticks
axes(a)          % set original axes as active
linkaxes([a b])  % link axes in case of zooming

%xlabel('Relative Days','FontWeight','b','FontSize',32)
%grid on
lblf2={'Watch Off';'Sleep';'Awake'};
%set(h35fpp,'YTick',.3:.7:2.3,'YTickLabel',lblf2,'FontSize',32)
%% Save figures
savefig(strcat(outp,stdy,'-',sb1(1:3),'-geneactiv_mtl3ss.fig'))
img = getframe(gcf);
imwrite(img.cdata, strcat(outp,stdy,'-',sb1(1:3),'-geneactiv_mtl3ss.png'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Keep this part  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except s

display('COMPLETE');
exit(0);
