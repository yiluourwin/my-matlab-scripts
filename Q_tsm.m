% ----------------------------------------------------------------
% Two Station Method for Q
% 
% Author: Yi LUO
% Built: 18-10-30 14:28
% Last Modified: 19-07-13 17:00
% ----------------------------------------------------------------
% Using Two Station Method to calculate Q. 
% Q=-pi.*frange.*(dist1-dist2).*u_avg.^-1.*(log((A1./A2).*(dist1./dist2).^0.5)).^-1;
% ----------------------------------------------------------------

%%
function [Q0,eta,param]=Q_tsm(sachd1,sacdata1,sachd2,sacdata2,frequency_limit,gv,noiseflag,figflag)


%% Input
% First four input arguements are required. 
% Records begin from event time. Using sachd.delta, sachd.b, sachd.dist. 
% Edit sachd.b to fix time. 
% The shorter record will be attached with 0. 
%% Default Input 
% frequency_limit=[0.1,2.5];    % limit output frequency. 
% gv=[3,3.5];   % gv is group velosity for target waveform, gv(1)<gv(2).
% noiseflag=0;  % if==1, reduce noise using waveform just before arriving time.
% figflag=0;    % if==1, print figures.
%% Default factors
% v_pn=8.5;     % pn phase velocity. Used in noise reduction. 
% df=0.02;      % moving smooth % log10(f)-log10(f1)=df
% dsmp_n=20;    % downsample factor. Only affects figures. 
% snr=2;        % Signal-noise rate. Signal lower this rate will be NaN.
% ----------------------------------------------------------------
%% Output
% param.frange      param.Q
% [raw amplitude]   [smoothed amp]      [result amp]
% param.A1          param.A_smooth1     param.A_result1   
% param.A2          param.A_smooth2     param.A_result2
%% if noiseflag
% [raw amplitude]   [smoothed amp]      [(snr<2)=NaN]   [length(snr>2)]
% param.NA1         param.NA_smooth1    param.NNA1      param.snlr1 
% param.NA2         param.NA_smooth2	param.NNA2      param.snlr2
% ----------------------------------------------------------------

%%
% defalt factors -----------------------------------
if nargin<5; frequency_limit=[0.1,2.5]; end
if nargin<6; gv=[3,3.5]; end % For Lg-wave
if nargin<7; noiseflag=0; end
if nargin<8; figflag=0; end

if isempty(frequency_limit); frequency_limit=[0.1,2.5]; end
if isempty(gv); gv=[3,3.5]; end % For Lg-wave
if isempty(noiseflag); noiseflag=0; end
if isempty(figflag); figflag=0; end
% --------------------------------------------------
v_pn=8.5;% pn phase velocity. Used in noise reduction. 
df=0.02; % moving smooth % log10(f)-log10(f1)=df % Reference 0.02
dsmp_n=20; % downsample factor. Only affects figures. 
snr=2; % Signal-noise rate. Signal lower this rate will be NaN. 
% --------------------------------------------------

warning('off','MATLAB:Axes:NegativeDataInLogAxis')

%%
if figflag
    figure(figflag)
end

% Data have different delta. Downsample the one with smaller delta.
if sachd1.delta>sachd2.delta
    sacdata2=resample(sacdata2,round(1/sachd1.delta),round(1/sachd2.delta));
    sachd2.delta=sachd1.delta;
elseif sachd2.delta>sachd1.delta
    sacdata1=resample(sacdata1,round(1/sachd2.delta),round(1/sachd1.delta));
    sachd1.delta=sachd2.delta;
end


% Data have different length. Add zeros to shorter one. 
dlength=length(sacdata1)-length(sacdata2);
if dlength>0
    sacdata2=[sacdata2;zeros(dlength,1)];
elseif dlength<0
    sacdata1=[sacdata1;zeros(-dlength,1)];
end

[param.A1,param.frange]=spectrum_tsm(sachd1,sacdata1,frequency_limit,gv,1,figflag);
[param.A2,~]=spectrum_tsm(sachd2,sacdata2,frequency_limit,gv,2,figflag);
    
param.A_smooth1=moving_smooth_tsm(param.A1,param.frange,df);
param.A_smooth2=moving_smooth_tsm(param.A2,param.frange,df);

if noiseflag
    [param.NA1,~]=noise_spectrum_tsm(sachd1,sacdata1,frequency_limit,gv,1,figflag,v_pn);
    [param.NA2,~]=noise_spectrum_tsm(sachd2,sacdata2,frequency_limit,gv,2,figflag,v_pn);
    param.NA_smooth1=moving_smooth_tsm(param.NA1,param.frange,df);
    param.NA_smooth2=moving_smooth_tsm(param.NA2,param.frange,df);
    
    param.NNA1=param.NA_smooth1;
    param.NNA1(param.A_smooth1./param.NA_smooth1<snr)=NaN;
    param.snlr1=sum(param.A_smooth1./param.NA_smooth1>snr);
    param.NNA2=param.NA_smooth2;
    param.NNA2(param.A_smooth2./param.NA_smooth2<snr)=NaN;
    param.snlr2=sum(param.A_smooth2./param.NA_smooth2>snr);
    
    param.A_result1=(param.A_smooth1.^2-param.NNA1.^2).^0.5;
    param.A_result2=(param.A_smooth2.^2-param.NNA2.^2).^0.5;
else
    param.A_result1=param.A_smooth1;
    param.A_result2=param.A_smooth2;
end

[param.Q,eta,Q0] = main_tsm (param.A_result1,sachd1.dist,param.A_result2,sachd2.dist,param.frange);


if figflag
    % 1
    sp_no(1)=subplot(2,3,1);
    title(sp_no(1),{'Record of St.1';['Distance=',num2str(sachd1.dist),'km']})
    xlabel('time / s')
    legend('show','autoupdate','off')
    for v_phases=[8.5,7,5,4,3.5,3,2.5,2]
        plot([sachd1.dist/v_phases,sachd1.dist/v_phases],[0,-max(abs(sacdata1))],'--','color',[0.466,0.674,0.188])
        text(sachd1.dist/v_phases,-max(abs(sacdata1)),num2str(v_phases))
    end
    xlim([sachd1.dist/10,sachd1.dist/1.2])
    hold off

    % 2
    sp_no(2)=subplot(2,3,2);
    title(sp_no(2),{'Record of St.2';['Distance=',num2str(sachd2.dist),'km']})
    xlabel('time / s')
    legend('show','autoupdate','off')
    for v_phases=[8.5,7,5,4,3.5,3,2.5,2]
        plot([sachd2.dist/v_phases,sachd2.dist/v_phases],[0,-max(abs(sacdata2))],'--','color',[0.466,0.674,0.188])
        text(sachd2.dist/v_phases,-max(abs(sacdata2)),num2str(v_phases))
    end
    xlim([sachd2.dist/10,sachd2.dist/1.2])
    hold off

	% 3
    sp_no(3)=subplot(2,3,3);
    loglog(param.frange,param.A1,':','color',[0,0.447,0.741],'displayname','St1 Raw')
    hold on
    loglog(param.frange,param.A_result1,'linewidth',2,'color',[0,0.447,0.741],'displayname','St1 Processed')
%     hold on
    loglog(param.frange,param.A2,':','color',[0.85,0.325,0.098],'displayname','St2 Raw')
    loglog(param.frange,param.A_result2,'linewidth',2,'color',[0.85,0.325,0.098],'displayname','St2 Processed')
    legend('show','location','best')
    xlabel('f / Hz');
    title(sp_no(3),'Spectra of 2 Records')
    hold off
    
    % 4
    sp_no(4)=subplot(2,3,4);
    loglog(param.frange,param.A_result1,'linewidth',2,'color',[0,0.447,0.741],'displayname','Processed Result')
    hold on
    loglog(param.frange,param.A1,'color',[0,0.447,0.741],'displayname','Raw Record')
    loglog(downsample(param.frange,dsmp_n),downsample(param.A_smooth1,dsmp_n),'o','markersize',4,'markerfacecolor','w','color',[0,0.447,0.741],'displayname','Moving Avg Record')
    if noiseflag
        loglog(param.frange,param.NA1,':','color',[0.929,0.694,0.125],'displayname','Noise') 
        loglog(downsample(param.frange,dsmp_n),downsample(param.NA_smooth1,dsmp_n),'^','markersize',4,'markerfacecolor','w','color',[0.929,0.694,0.125],'displayname','Moving Avg Noise')
    end
    legend('show','location','best')
    xlabel('f / Hz');
    title(sp_no(4),'Spectra of St.1')
    hold off
    
    % 5
    sp_no(5)=subplot(2,3,5);
    loglog(param.frange,param.A_result2,'linewidth',2,'color',[0.85,0.325,0.098],'displayname','Processed Result')
    hold on
    loglog(param.frange,param.A2,'color',[0.85,0.325,0.098],'displayname','Raw Record')
    loglog(downsample(param.frange,dsmp_n),downsample(param.A_smooth2,dsmp_n),'o','markersize',4,'markerfacecolor','w','color',[0.85,0.325,0.098],'displayname','Moving Avg Record')
    if noiseflag
        loglog(param.frange,param.NA2,':','color',[0.929,0.694,0.125],'displayname','Noise')
        loglog(downsample(param.frange,dsmp_n),downsample(param.NA_smooth2,dsmp_n),'^','markersize',4,'markerfacecolor','w','color',[0.929,0.694,0.125],'displayname','Moving Avg Noise')
    end
    legend('show','location','best')
    xlabel('f / Hz');
    title(sp_no(5),'Spectra of St.2')
    hold off
    
    % 6
    sp_no(6)=subplot(2,3,6);
    % yyaxis left
    % set(gca,'ycolor','k','XScale','log','YScale','log')
    
    loglog(param.frange,param.Q);
    
    hold on
    % yyaxis right
    
    loglog(param.frange,Q0.*param.frange.^eta);
    
    % ylim([100,2000]);
    xlim([param.frange(1),param.frange(end)]);
    grid on
    
    % set(gca,'xtick',[0.25 0.5 1 2],'xticklabel',[0.25 0.5 1 2],...
    %     'ytick',[100 250 500 750 1000 1500],'yticklabel',[100 250 500 750 1000 1500]);
    
    xlabel('f / Hz');
    ylabel('Q');
    hold off
    
    title(sp_no(6),{'Q - f';['Q_0=',num2str(Q0),'    \eta=',num2str(eta)]});
end

% param.frange=frange;
% param.Q=Q;
% 
% param.A1=A1;
% %param.frange1=frange1;
% param.A_smooth1=A_smooth1;
% param.A_result1=A_result1;
% 
% param.A2=A2;
% %param.frange2=frange2;
% param.A_smooth2=A_smooth2;
% param.A_result2=A_result2;
% 
% if noiseflag
%     param.NA1=NA1;
%     param.NA_smooth1=NA_smooth1;
%     param.NA2=NA2;
%     param.NA_smooth2=NA_smooth2;
% end

end

%%
% --------------------------------
% Calculate spectrum
% Two Station Method
% sacfile begins from event time
% Built: 18-10-29 14:18
% Last Modified: 18-11-01 21:24
% --------------------------------

function [amplitude,frange]=spectrum_tsm(sachd,sacdata,frequency_limit,gv,fig_no,debug)

% fmax=2.5;
% fmin=0.1;

% gv=[3,3.5];
% R=6371;

% [sachd,sacdata]=rsac(sacfile);

Fs=1./sachd.delta;
%dist=sachd.dist;

dt=2;% Window slip
tt=1:size(sacdata,1);
tt=heaviside(tt-(sachd.dist./gv(2)-sachd.b).*Fs).*heaviside((sachd.dist./gv(1)-sachd.b).*Fs-tt)...
    +heaviside(tt-(sachd.dist./gv(2)-sachd.b-dt).*Fs).*heaviside((sachd.dist./gv(2)-sachd.b).*Fs-tt)...
    .*(1./(dt.*Fs).*(tt-(sachd.dist./gv(2)-sachd.b-dt).*Fs))...
    +heaviside(tt-(sachd.dist./gv(1)-sachd.b).*Fs).*heaviside((sachd.dist./gv(1)-sachd.b+dt).*Fs-tt)...
    .*(-1./(dt.*Fs).*(tt-(sachd.dist./gv(1)-sachd.b+dt).*Fs));% Assume that all the sachd.kztime equal. Actrually some differs about 10ms.
% tt=heaviside(tt-(dist./gv(2)-t(1)).*Fs).*heaviside((dist./gv(1)-t(1)).*Fs-tt)...
%     +heaviside(tt-(dist./gv(2)-t(1)-dt).*Fs).*heaviside((dist./gv(2)-t(1)).*Fs-tt)...
%     .*(1./(dt.*Fs).*(tt-(dist./gv(2)-t(1)-dt).*Fs))...
%     +heaviside(tt-(dist./gv(1)-t(1)).*Fs).*heaviside((dist./gv(1)-t(1)+dt).*Fs-tt)...
%     .*(-1./(dt.*Fs).*(tt-(dist./gv(1)-t(1)+dt).*Fs));
tw=diag(tt);
data=tw*sacdata;

n=size(data,1);

n=2^nextpow2(n);

fdata=fft(data,n,1);

P2=abs(fdata./n);
P1=P2(1:n/2+1);
P1(2:end-1)=2*P1(2:end-1);

frange=Fs*(0:(n/2))/n;
frange=frange';
amplitude=P1;

amplitude=amplitude(frange>=frequency_limit(1) & frange<=frequency_limit(2));
frange=frange(frange>=frequency_limit(1) & frange<=frequency_limit(2));

if debug
    subplot(2,3,fig_no)
    plot((0:size(sacdata)-1).*sachd.delta+sachd.b,sacdata,'displayname','Record')
    hold on
    plot((0:size(sacdata)-1).*sachd.delta+sachd.b,data,'displayname','Selected Phase')
%     plot([0,dist./gv(2)-dt,dist./gv(2),dist./gv(1),dist./gv(1)+dt,size(sacdata,1).*sachd.delta],...
%         [0,0,max(sacdata),max(sacdata),0,0],'r');
    
end

end

%%
% --------------------------------
% Calculate noise spectrum
% Two Station Method
% sacfile begins from event time
% Built: 18-11-05 15:30
% Last Modified: 18-11-05 15:43
% --------------------------------

function [noise_amplitude,frange]=noise_spectrum_tsm(sachd,sacdata,frequency_limit,gv,fig_no,debug,v_pn)

% fmax=2.5;
% fmin=0.1;

% gv=[3,3.5];
% R=6371;

%v_pn=8.5;

% [sachd,sacdata]=rsac(sacfile);

Fs=1./sachd.delta;
%dist=sachd.dist;

t_move=-sachd.dist./gv(1)+sachd.dist./v_pn;

dt=2;% Window slip
tt=1:size(sacdata,1);
tt=heaviside(tt-(sachd.dist./gv(2)-sachd.b+t_move).*Fs).*heaviside((sachd.dist./gv(1)-sachd.b+t_move).*Fs-tt)...
    +heaviside(tt-(sachd.dist./gv(2)-sachd.b-dt+t_move).*Fs).*heaviside((sachd.dist./gv(2)-sachd.b+t_move).*Fs-tt)...
    .*(1./(dt.*Fs).*(tt-(sachd.dist./gv(2)-sachd.b-dt+t_move).*Fs))...
    +heaviside(tt-(sachd.dist./gv(1)-sachd.b+t_move).*Fs).*heaviside((sachd.dist./gv(1)-sachd.b+dt+t_move).*Fs-tt)...
    .*(-1./(dt.*Fs).*(tt-(sachd.dist./gv(1)-sachd.b+dt+t_move).*Fs));% Assume that all the sachd.kztime equal. Actrually some differs about 10ms.
% tt=heaviside(tt-(dist./gv(2)-t(1)).*Fs).*heaviside((dist./gv(1)-t(1)).*Fs-tt)...
%     +heaviside(tt-(dist./gv(2)-t(1)-dt).*Fs).*heaviside((dist./gv(2)-t(1)).*Fs-tt)...
%     .*(1./(dt.*Fs).*(tt-(dist./gv(2)-t(1)-dt).*Fs))...
%     +heaviside(tt-(dist./gv(1)-t(1)).*Fs).*heaviside((dist./gv(1)-t(1)+dt).*Fs-tt)...
%     .*(-1./(dt.*Fs).*(tt-(dist./gv(1)-t(1)+dt).*Fs));
tw=diag(tt);
data=tw*sacdata;

n=size(data,1);

n=2^nextpow2(n);

fdata=fft(data,n,1);

P2=abs(fdata./n);
P1=P2(1:n/2+1);
P1(2:end-1)=2*P1(2:end-1);

frange=Fs*(0:(n/2))/n;
frange=frange';
noise_amplitude=P1;

% noise_amplitude=noise_amplitude(frange>=fmin & frange<=fmax,:);
% frange=frange(frange>=fmin & frange<=fmax);
noise_amplitude=noise_amplitude(frange>=frequency_limit(1) & frange<=frequency_limit(2));
frange=frange(frange>=frequency_limit(1) & frange<=frequency_limit(2));

if debug
    subplot(2,3,fig_no)
%     legend('autoupdate','on')
    plot((0:size(sacdata)-1).*sachd.delta+sachd.b,data,'color',[0.929,0.694,0.125],'displayname','Noise')
%     plot([0,dist./gv(2)-dt+t_move,dist./gv(2)+t_move,dist./gv(1)+t_move,dist./gv(1)+dt+t_move,size(sacdata,1).*sachd.delta],...
%         [0,0,max(sacdata),max(sacdata),0,0],'y');
%     hold off
end

end

%%
% --------------------------------
% Main part
% Two Station Method
%
% Built: 18-10-30 14:14
% Last Modified: 18-10-30 14:28
% --------------------------------

function [Q,eta,Q0] = main_tsm (A1,dist1,A2,dist2,frange)

u_avg=3.3;

Q=-pi.*frange.*(dist1-dist2).*u_avg.^-1.*(log((A1./A2).*(dist1./dist2).^0.5)).^-1;

p=polyfit(log10(frange(Q>0)),log10(Q(Q>0)),1);
Q0=10.^p(2);
eta=p(1);

end

%%
% --------------------------------
% moving smooth: root mean square
% Two Station Method
%
% Built: 18-12-17 14:14
% Last Modified: 18-12-17 14:28
% --------------------------------

function [ax]=moving_smooth_tsm(a1,f,df)

ax=zeros(length(a1),1);
for i=1:length(f)
    ax(i)=sqrt(mean((a1(f>(f(i).*10.^-df) & f<(f(i).*10.^+df))).^2));
    % log10(f)-log10(f1)=df
end

%a=interp1(f1,ax,f);

end