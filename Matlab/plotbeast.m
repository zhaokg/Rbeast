function H = plotbeast(o, varargin )
%   USAGE: <strong>plotbeast(o, ...) </strong>
%
%   <strong>o</strong>:  the time series analysis output from  beast123; o
%   should contain results for mulltiple time series
%
%   <strong>'index'</strong>: if o contains results for more than 1 time
%   series, index specifies for which time series the result is printed.
%   If o is the result for a 3D stacked cube, index will be a vector of 2
%   integer to specify the row and col of the desired pixel. If o contains
%   only one time series, index will be ignored
%
%   <strong>'vars'</strong>:a vector of strings indicating the elements or variables 
%   of o to plot. Possible vars strings include 'st' (season plus trend), 's' (season component),
%   't' (trend component), 'o' (outliers), 'scp', 'tcp', 'ocp'  (occurrence probability
%   of seasonal/trend/outlier changepoint), 'sorder' (seasonal harmonic order), 
%   'torder' (trend polynomial order), 'samp' (amplitude of seasonality), 'tslp' (slope of trend),
%   'slpsgn' (probabilities of the slope being positive, zero, and negative) and 'error' (remainder).
%  
%   <strong>'ncpStat'</strong>: A string to specify which statistic is used for 
%   the Number of ChangePoint (ncp). Five values are possible: 'mean', 'mode', 'median',
%   'pct90', and 'max'; the default is 'median'. Individual models sampled by BEAST has
%   a varying dimension (e.g., number of changepoints or knots). For example, 
%   if mcmc$samples=10, the numbers of changepoints for the 10 sampled models are assumed
%   to be c(0, 2, 4, 1, 1, 2, 7, 6, 6, 1). The mean ncp will be 3.1 (rounded to 3), 
%   the median is 2.5 (2), the mode is 1, and the maximum is 7.  The 'max' option plots 
%   all the changepoints recorded in out$trend$cp, out$season$cp, or out$outlier$cp; many of 
%   these changepoints are bound to be false positives, so do not blindly treat all of them 
%   as actual changepoints.
%
%   <strong>Examples</strong>:
%   % in the slpsgn plot, the upper red envelope refers to probability of trend slope being positive,
%   % the middle green to probability of trend slope being zero, the blue to probability of slope
%   % of being negative
%   plotbeast( beast(beach), 'vars',["st", "t", "tcp" ,"slpsgn"], "ncpStat",'mean')
% 
%   plotbeast(o, 'index',3)
%   plotbeast(o, 'index',[3,5])
%
%   <strong>Contact info</strong>: To report bug or get help, do not hesitate to contact Kaiguang Zhao
%   at <strong>zhao.1423@osu.edu</strong>.
%
%%
if ~strcmp(o.class,'beast')
    error('the input has to be an output from the BEAST functions');
end

%% Check the second argument -- the option parameter
    n=length(varargin);
    if mod(n,2)~=0
        msg=[ "the optional arg list must be paired keywords and values; the number of extra args must be even.\n", ...
            " Examples:  plotbeast(y) \n", ...
            "            plotbeast(y,'vars',""st"",""s""],'ncpStat','median') \n", ...
            ];
        msg=sprintf(strcat(msg{:}));
        error(msg);
    end

    KeyList   = varargin(1:2:n);
    KeyList   = cellfun(@char,KeyList,'UniformOutput',false); % convert strings--if any--to chars,
    ValList   = varargin(2:2:n);

    indexDefautValue    = 1;
    varsDefautValue     = ["st","s","scp","sorder","t","tcp","torder", "slpsgn","o","ocp","error"];
    ncpStatDefaultValue = 'median';
    index   = GetValueByKey(KeyList, ValList, 'index',   indexDefautValue);
    vars    = GetValueByKey(KeyList, ValList, 'vars',    varsDefautValue);
    ncpStat = GetValueByKey(KeyList, ValList, 'ncpStat', ncpStatDefaultValue);
    opt     = GetValueByKey(KeyList, ValList, 'layout',  []);
    if isempty(opt)
        opt.leftmargin    =0.1;
        opt.rightmargin   =0.1;
        opt.topmargin     =0.05;
        opt.bottommargin  =.07;
        opt.verticalspace =0.01;
    end
    fig  = GetValueByKey(KeyList, ValList, 'fig', []);
%%
vars     = lower(vars);
vars_log = vars=="st"|vars=="s"|vars=="t"|vars=="scp"|vars=="tcp" ...
    |vars=="sorder"|vars=="torder"|vars=="error"|vars=="o"|vars=="ocp"|vars=="samp"|vars=="tslp"|vars=="slpsgn";
vars     = vars(vars_log);

%%
ylab     = vars;
ylab(vars=='st')  ='Y';
ylab(vars=='s')   ='season';
ylab(vars=='t')   ='trend';
ylab(vars=='o')   ='outlier';
ylab(vars=='scp') ="Pr(scp)";
ylab(vars=='tcp') ="Pr(tcp)";
ylab(vars=='ocp') ="Pr(ocp)";
ylab(vars=='sorder')= "sOrder";
ylab(vars=='torder')= "tOrder";
ylab(vars=='samp')  = 'amplitude';
ylab(vars=='tslp')  = 'slope';
ylab(vars=='slpsgn') = "slpSign";
ylab(vars=='error')  = "error";

col=cell(1,length(vars));
col(vars=='st')    ={[0.1,0.1,0.1]};
col(vars=='s')     ={'r'};
col(vars=='scp')   ={'r'};
col(vars=='sorder')={'r'};
col(vars=='samp')  ={'r'};
col(vars=='t')     ={'g'};
col(vars=='tcp')   ={'g'};
col(vars=='torder')={'g'};
col(vars=='tslp')  ={'g'};
col(vars=='slpsgn')={'k'} ; %the border only
col(vars=='o')     ={'b'};
col(vars=='ocp')   ={'b'};
col(vars=='error') ={[0.4,0.4,0.4]};

heights=zeros(1,length(vars))+1;
heights(vars=='st')=.8;
heights(vars=='s')= .8;
heights(vars=='t')=.8;
heights(vars=='o')=.8;
heights(vars=='scp')=.5;
heights(vars=='tcp')=.5;
heights(vars=='ocp')=.5;
heights(vars=='sorder')=.5;
heights(vars=='torder')=.5;
heights(vars=='samp')   =.4;
heights(vars=='tslp')   =.4;
heights(vars=='slpsgn') =.4;
heights(vars=='error')=.4;
%%
if ( length(o.marg_lik)> 1 )
    %# more than time series is present
    o=extractbeast(o,index);
end
if (isnan(o.marg_lik))
    H=[];
    warning("The result of the selected time series selected is invalid.");
    return;
end
%%
hasData    = isfield(o,'data') && ~isempty(getfield(o,'data'));
hasSeason  = isfield(o,'season') &&  ~isempty(o.season);
if hasSeason 
    hasSOrder = isfield(o.season,'order') &&  ~isempty(o.season.order);
    hasAmp    = isfield(o.season,'amp')   &&  ~isempty(o.season.amp);
else
    hasSOrder = false;
    hasAmp    = false;
end

hasOutlier  = isfield(o,'outlier')    && ~isempty(o.outlier);
hasTOrder   = isfield(o.trend,'order') &&  ~isempty(o.trend.order);
hasSlp      = isfield(o.trend,'slp')   &&  ~isempty(o.trend.slp);
 
has=[];
has.hasAmp=hasAmp;
has.hasSlp=hasSlp;
has.hasSeason=hasSeason;
has.hasSOrder=hasSOrder;
has.hasTOrder=hasTOrder;
has.hasOutlier=hasOutlier;
has.hasData=hasData;
%%
idx = (1:length(vars)) > 0;
if(~hasAmp)       idx   = idx & ~(vars=='samp');         end
if(~hasSlp)       idx   = idx & ~(vars=='tslp'|vars=='slpsgn'); end
if(~hasSeason)    idx   = idx & ~(vars=='st'|vars=='s'|vars=='sorder'|vars=='scp') ; end
if(~hasSOrder)    idx   = idx & ~(vars=='sorder') ; end
if(~hasTOrder)    idx   = idx & ~(vars=='torder') ; end
if(~hasOutlier)   idx   = idx &~(vars=='o'|vars=='ocp'); end
if(~hasData)      idx   = idx &~(vars=='error') ;  end

col         = col(idx);
vars        = vars(idx);
ylab        = ylab(idx);
heights     = heights(idx);

nPlots=length(vars);
if(nPlots==0)
    error("No valid variable names speciffied int the 'vars' argument. Possible names include 'st','t','s','sorder','torder','scp','tcp','samp','tslp','o', 'ocp', and 'error'. ");
end
%%


%#######################################################
%#  Functions and variables to load the outputs
%########################################################

if isempty(fig)    
    clf;
    fig=gcf;
end
H= axeslayout(opt, heights,fig);

% #######################################################
% #  Create a subplot given the relative heights of the vertical plots
% ########################################################
x     = o;

t     = x.time;
t2t   = [t; t(end:-1:1)];
N     = length(t);
for i =1:length(vars)
    
    ytitle = ylab(i);
    var    = vars(i);
    clr    = col{i};
    
    h = H(i);
    axes(h);    
    cla(h);         % cla
    hold(h,'on');   % hold on
 
    if (var=='st')
        [Yts,YtsSD, Yerr]=  get_Yts(x,hasSeason,hasOutlier, hasData);
        plot_st(h,ytitle, has,clr,x,t,t2t,Yts,YtsSD) ;
    end
    if (var=='s' )
        [Y,SD,CI, Amp,AmpSD, Order] = get_S(x,hasAmp,hasSOrder);
        [cp, cpCI, ncp, ncpPr,cpPr, cpChange, Prob, Prob1] = get_scp(x,ncpStat);
        plot_y(h,ytitle, has,clr,x, t,t2t,Y, CI, ncp, cp) ;
    end
    if (var=='t' )
        [Y,SD, CI,Slp,SlpSD,SlpSignPos,SlpSignZero,Order] = get_T(x,hasSlp,hasTOrder);
        [cp, cpCI, ncp, ncpPr,cpPr, cpChange, Prob, Prob1] = get_tcp(x,ncpStat);
        
        plot_y(h,ytitle, has,clr,x, t,t2t,Y, CI, ncp, cp) ;
    end
    
    if (var=='scp' )
        [Y,SD,CI, Amp,AmpSD, Order] = get_S(x,hasAmp,hasSOrder);
        [cp, cpCI, ncp, ncpPr,cpPr, cpChange, Prob, Prob1] = get_scp(x,ncpStat);
        plot_prob(h,ytitle,  has, clr,x, t, t2t, Prob1, Prob,ncp,cp );
    end
    
    if (var=='tcp' )
        [Y,SD, CI,Slp,SlpSD,SlpSignPos,SlpSignZero,Order] = get_T(x,hasSlp,hasTOrder);
        [cp, cpCI, ncp, ncpPr,cpPr, cpChange, Prob, Prob1] = get_tcp(x,ncpStat);
        plot_prob(h,ytitle,  has, clr,x, t, t2t, Prob1, Prob,ncp,cp );
        
    end

    if (var=='sorder' )
        [Y,SD,CI, Amp,AmpSD, Order] = get_S(x,hasAmp,hasSOrder);
        [cp, cpCI, ncp, ncpPr,cpPr, cpChange, Prob, Prob1] = get_scp(x,ncpStat);
       plot_order( h,ytitle, has, clr,x, t, t2t,Order, ncp, cp);
        
    end
    
    if (var=='torder' )
        [Y,SD, CI,Slp,SlpSD,SlpSignPos,SlpSignZero,Order] = get_T(x,hasSlp,hasTOrder);
        [cp, cpCI, ncp, ncpPr,cpPr, cpChange, Prob, Prob1] = get_tcp(x,ncpStat);
       plot_order(h,ytitle, has, clr,x, t, t2t,Order, ncp, cp);
        
    end    
    
    if (var=='samp' )
        [Y,SD,CI, Amp,AmpSD, Order] = get_S(x,hasAmp,hasSOrder);
        plot_amp(h,ytitle,has, clr,x, t, t2t,Amp, AmpSD);
        
    end
    
    if (var=='tslp' )
        [Y,SD, CI,Slp,SlpSD,SlpSignPos,SlpSignZero,Order] = get_T(x,hasSlp,hasTOrder);
        plot_slp( h,ytitle,has, clr,x, t, t2t,Slp,SlpSD)
        
    end
    
    if (var=='slpsgn' )
        [Y,SD, CI,Slp,SlpSD,SlpSignPos,SlpSignZero,Order] = get_T(x,hasSlp,hasTOrder);
        plot_slpsgn(h,ytitle, has, clr,x, t, t2t,SlpSignPos,SlpSignZero)
        
    end
   
    if (var=='o' )
        [Y,SD,CI ]  = get_O(x );
        plot_o( h,ytitle,has, clr,x, t, t2t,Y,SD);
        
    end
    
    if (var=='ocp' )
       [cp, cpCI, ncp, ncpPr,cpPr, cpChange, Prob, Prob1] = get_ocp(x,ncpStat);
       plot_oprob(h,ytitle, has, clr,x, t, t2t, Prob1, Prob,ncp,cp );        
    end
    
    if (var=='error' )        
        [Yts,YtsSD, Yerr ] = get_Yts(x,hasSeason,hasOutlier, hasData);        
        plot_error(h,ytitle, has, clr,x, t, t2t, Yerr);        
    end   
    
    if mod(i,2)==1
        set(h,'YaxisLocation','left');
    else
        set(h,'YaxisLocation','right');
    end
    
    set(h,'YticklabelRotation',90);    
      
    if(i==1)
     title(h,'BEAST decompositon and changepoint detection');
    end

    if(i==nPlots)
	  xlabel(h,'time');
    else
      set(h,'xticklabel',[]);
    end   
   
    ylabel(h,ytitle);
    
end
end
%%

%% Functions to return a default value if the field is missing from opt
function value=GetValueByKey(KeyList, ValList, key,defaultValue)
   idx=find(strcmp(KeyList,key));
   if isempty(idx)
       value=defaultValue;
   else
       value=ValList{idx(1)};
   end
end
%%
function H = axeslayout(opt, hLayout, fig)

lm=opt.leftmargin;
rm=opt.rightmargin;
tm=opt.topmargin;
bm=opt.bottommargin;
vd=opt.verticalspace;

%hLayout=[1 1 -1 1 1 -1 1 1];
winList=find(hLayout>0);
hgt=(1 - tm -bm -abs( sum(hLayout(hLayout<0))*vd ))/sum(hLayout(winList));
H  =[];
for i=1:length(winList)
    H(i) = axes(fig);
    idx  =  winList(i);
    if idx==1
        yup=1-tm;
    else
        slist=hLayout(1:idx-1);
        slist(slist >0)=slist(slist >0)*hgt;
        slist(slist <0)=-slist(slist <0)*vd;
        yup=1-tm-sum(slist);
    end
    
    set(H(i),'pos',[lm, yup-hLayout(idx)*hgt, 1-rm-lm,hLayout(idx)*hgt ]);
    set(H(i),'box','on');
end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [Yts,YtsSD, Yerr ] = get_Yts(x,hasSeason,hasOutlier, hasData)

Yts   = x.trend.Y;
SD2   = x.trend.SD.^2 +  x.sig2(1);

if hasSeason
    Yts  =  Yts+x.season.Y;
    SD2  =  SD2+x.season.SD.^2;
end

if hasOutlier
    Yts  =  Yts+x.outlier.Y;
    SD2  =  SD2+x.outlier.SD.^2;
end

SD    = sqrt(SD2);
tmp   = Yts+SD;
YtsSD = [Yts-SD; tmp(end:-1:1)];

if hasData
    Yerr  = x.data-Yts;
else
    Yerr=[];
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [Y,SD, CI,Slp,SlpSD,SlpSignPos,SlpSignZero,Order] = get_T(x,hasSlp,hasTOrder)

Y  = x.trend.Y;
tmp=Y+x.trend.SD;
SD =[Y-x.trend.SD;  tmp(end:-1:1)];
if isfield(x.trend,'CI') &  ~isempty(x.trend.CI)
    tmp = x.trend.CI(:,2);
    CI  = [x.trend.CI(:,1); tmp(end:-1:1)] ;
else
    CI  =SD;
end

if (hasSlp)
    Slp         =  x.trend.slp;
    tmp         =  Slp+x.trend.slpSD;
    SlpSD       =  [Slp-x.trend.slpSD;  tmp(end:-1:1)];
    SlpSignPos  =  x.trend.slpSgnPosPr;
    SlpSignZero =  x.trend.slpSgnZeroPr;
else
    Slp=[];
    SlpSD=[];
    SlpSignPos=[];
    SlpSignZero=[];
end

if (hasTOrder)
    Order       =  x.trend.order;
else
    Order=[];
end


end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [cp, cpCI, ncp, ncpPr,cpPr, cpChange, Prob, Prob1] = get_tcp(x,ncpStat)

cmpnt   =  x.trend;

cp     = cmpnt.cp;
cpCI   = cmpnt.cpCI;
 
switch lower(ncpStat)
    case 'mode'
        ncp=cmpnt.ncp_mode;
    case 'median'
        ncp= cmpnt.ncp_median;
    case 'mean'
        ncp=cmpnt.ncp;
    case 'pct90'
        ncp=cmpnt.ncp_pct90;
    case 'pct10'
        ncp=cmpnt.ncp_pct10;
    case 'max'
        ncp=sum(~isnan(cp));
    otherwise
        ncp=cmpnt.ncp_mode;
end

ncp    = round(ncp);
ncpPr  = cmpnt.ncpPr;
cpPr   = cmpnt.cpPr;
cpChange= cmpnt.cpAbruptChange;

Prob    = cmpnt.cpOccPr;
Prob1   = [Prob;Prob-Prob];
end
%% ###########################################################
function  [Y,SD,CI, Amp,AmpSD, Order] = get_S(x,hasAmp,hasSOrder)
Y      = x.season.Y;
tmp    = Y+x.season.SD;
SD     =[Y-x.season.SD;  tmp(end:-1:1)];

if isfield(x.season,'CI') &  ~isempty(x.season.CI)
    tmp = x.season.CI(:,2);
    CI  = [x.season.CI(:,1); tmp(end:-1:1)] ;
else
    CI  =SD;
end

if hasAmp
    Amp    =x.season.amp;
    tmp    =Amp+x.season.ampSD;
    AmpSD  =[Amp-x.season.ampSD;tmp(end:-1:1) ];
else
    Amp=[];
    AmpSD=[];
end

if hasSOrder
    Order=x.season.order;
else
    Order=[];
end
end


%% ###########################################################
function  [Y,SD,CI ] = get_O(x)
Y      = x.outlier.Y;
tmp    = Y+x.outlier.SD;
SD     =[Y-x.outlier.SD;  tmp(end:-1:1)];

if isfield(x.outlier,'CI') &  ~isempty(x.outlier.CI)
    tmp = x.outlier.CI(:,2);
    CI  = [x.outlier.CI(:,1); tmp(end:-1:1)] ;
else
    CI  =SD;
end
end
 
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [cp, cpCI, ncp, ncpPr,cpPr, cpChange, Prob, Prob1] = get_scp(x,ncpStat)

cmpnt   =  x.season;

cp     = cmpnt.cp;
cpCI   = cmpnt.cpCI;

switch lower(ncpStat)
    case 'mode'
        ncp=cmpnt.ncp_mode;
    case 'median'
        ncp= cmpnt.ncp_median;
    case 'mean'
        ncp=cmpnt.ncp;
    case 'pct90'
        ncp=cmpnt.ncp_pct90;
    case 'pct10'
        ncp=cmpnt.ncp_pct10;
    case 'max'
        ncp=sum(~isnan(cp));
    otherwise
        ncp=cmpnt.ncp_mode;
end

ncp    = round(ncp);
ncpPr  = cmpnt.ncpPr;
cpPr   = cmpnt.cpPr;
cpChange= cmpnt.cpAbruptChange;

Prob    = cmpnt.cpOccPr;
Prob1   = [Prob;Prob-Prob];
end

function  [cp, cpCI, ncp, ncpPr,cpPr, cpChange, Prob, Prob1] = get_ocp(x,ncpStat)

cmpnt   =  x.outlier;

cp     = cmpnt.cp;
cpCI   = cmpnt.cpCI;

switch lower(ncpStat)
    case 'mode'
        ncp=cmpnt.ncp_mode;
    case 'median'
        ncp= cmpnt.ncp_median;
    case 'mean'
        ncp=cmpnt.ncp;
    case 'pct90'
        ncp=cmpnt.ncp_pct90;
    case 'pct10'
        ncp=cmpnt.ncp_pct10;
    case 'max'
        ncp=sum(~isnan(cp));
    otherwise
        ncp=cmpnt.ncp_mode;
end

ncp    = round(ncp);
ncpPr  = cmpnt.ncpPr;
cpPr   = cmpnt.cpPr;
cpChange= []; % no delta_change for the outlier component

Prob    = cmpnt.cpOccPr;
Prob1   = [Prob;Prob-Prob];
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_st(h,ytitle,has,clr,x, t,t2t,Yts,YtsSD)
 
alpha=0.2;
fill(h, t2t, YtsSD, clr,'LineStyle','none' ,'FaceAlpha',alpha);
if (has.hasData)
    plot(h,t, x.data,'o','color',clr);
end
plot(h, t, Yts,  'color',clr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_y(h,ytitle,has,clr,x, t,t2t,Y, CI, ncp, cp)

alpha=0.2;
if (has.hasData && ~has.hasSeason)
    plot(h,t, x.data,'o','color',[.5,.5,.5]);
end
fill(h, t2t, CI,clr,'LineStyle','none','FaceAlpha',0.1) ;

plot(h,t,Y,clr);
ylim=get(h,'ylim');
for i = 1 : ncp
    plot(h, [cp(i),cp(i)], ylim, 'color', 'k');
end

set(h,'ylim',ylim);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  plot_prob(h,ytitle, has, clr,x, t, t2t, Prob1, Prob,ncp,cp )

%cla;  hold on;

alpha=0.2;

fill(h, t2t, Prob1,  clr,'LineStyle','none','FaceAlpha',alpha) ;
plot(h,t, Prob, clr);
maxp=min(1,max(Prob)*1.5);
maxp=max(maxp,0.2);
set(h,'ylim',[0,maxp]);

for i = 1 : ncp
    plot( h,[cp(i),cp(i)], get(h,'Ylim'), 'color', 'k');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_order(h,ytitle, has, clr,x, t, t2t,Order, ncp, cp)
%cla;  hold on;
plot(h, t,   Order,  'color',clr );
maxp=max(max(Order),1.05);
minp=-0.05;
set(h,'ylim',[minp,maxp]);

for i = 1 : ncp
    plot(h, [cp(i),cp(i)], get(h,'Ylim'), 'color', 'k');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   plot_amp(h,ytitle,has, clr,x, t, t2t,Amp, AmpSD)
%cla;  hold on;
alpha=0.5 ;
%$fill(t2t, AmpSD,   col  = rgb(col[1],col[2],col[3],alpha), border = NA);
plot(h, t,   Amp,'color',clr );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_slp(h,ytitle, has, clr,x, t, t2t,Slp,SlpSD)
%cla;  hold on;
alpha=0.5;
fill(h,t2t, SlpSD,  clr,'LineStyle','none','FaceAlpha',alpha) ;
plot(h,t,Slp, 'color',clr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_slpsgn (h,ytitle,has, clr,x, t, t2t,SlpSignPos,SlpSignZero )
%cla;  hold on;
alpha=0.5;
SlpSignNeg = 1-SlpSignPos-SlpSignZero;

y2y   = [t-t; SlpSignNeg(end:-1:1)];           fill(t2t, y2y, [0,0,1], 'LineStyle','none','FaceAlpha',alpha) ;
y2y   = [SlpSignNeg;  1-SlpSignPos(end:-1:1)]; fill(t2t, y2y, [0,1,0], 'LineStyle','none','FaceAlpha',alpha) ;
y2y   =[ 1-SlpSignPos; t-t+1];                fill(t2t, y2y, [1,0,0], 'LineStyle','none','FaceAlpha',alpha) ;
plot(h, t,  t-t+0.5 );
end

function  plot_o(h,ytitle,has,clr,x, t,t2t,Y, CI, ncp, cp)
%cla;  hold on;
alpha=0.5;

%#polygon(t2t, CI,   col  = rgb(col[1],col[2],col[3],alpha), border = NA);
%#points( t,   Y,    type = 'l',col='#333333');
stem(h,t,Y,'marker','none','color',clr);
end

function   plot_oprob(h,ytitle, has, clr,x, t, t2t, Prob1, Prob,ncp,cp )
alpha=0.2
%plot( c(t2t[1],t2t), c(0.22,Prob1),type = 'n', ann=FALSE, xaxt='n', yaxt='n');
%# polygon(t2t, Prob1, col  = rgb(col[1],col[2],col[3],alpha), border = NA);
%# points( t,   Prob,  col  = rgb(col[1],col[2],col[3])  ,       lwd = 1,type = 'l' );
stem(h,t,Prob,'marker','none','color',clr);
end

function plot_error(h,ytitle,has, clr,x, t, t2t, Yerr)
%plot(   t,   Yerr,  type = 'n',ann = FALSE, xaxt = 'n', yaxt = 'n');
plot(h,t, t-t, 'color',clr);
stem(h,t,Yerr,'color',clr,'marker','none');
end
