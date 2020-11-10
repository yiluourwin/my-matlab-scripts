% --------------------------------
% calculate and plot spectra
% 
% Author: Yi LUO
% Built: 2019/06/01 16:00
% Last Modified: 2020/11/07 14:00
% --------------------------------

function [a,f,x]=spectra_cal(Fs,x,frange,df,twindow,figflag,plottype)
% use [] to ignore input arguements
% x is n-by-1(or m) metrix
% frange and twindow are 1-by-2 metrix. twindow relatives to the beginning 0 time. 
% If x has more than 1 column, they should share the same sample rate, frange and twindow. 
% Smoothing param df=0.2 is recommanded.
% plottype	defalt:loglog   options: plot, semilog


if nargin<3; frange=nan; end
if nargin<4; df=nan; end
if nargin<5; twindow=nan; end
if nargin<6; figflag=nan; end
if nargin<7; plottype='loglog';end

if isempty(frange); frange=nan; end
if isempty(df); df=nan; end
if isempty(twindow); twindow=nan; end
if isempty(figflag) || figflag==0; figflag=nan; end

%Fs=mean(1./mean(diff(t),1));

if ~isnan(twindow)
    st=(twindow(2)-twindow(1))/20; % Window slip second
    tt=1:size(x,1);
    tt=heaviside(tt-twindow(1).*Fs).*heaviside(twindow(2).*Fs-tt)...
        +heaviside(tt-(twindow(1)-st).*Fs).*heaviside(twindow(1).*Fs-tt)...
        .*(1./(st.*Fs).*(tt-(twindow(1)-st).*Fs))...
        +heaviside(tt-twindow(2).*Fs).*heaviside((twindow(2)+st).*Fs-tt)...
        .*(-1./(st.*Fs).*(tt-(twindow(2)+st).*Fs));
    tw=diag(tt);
    x=tw*x;
    
    L=Fs*(twindow(2)-twindow(1));
else
    L=size(x,1);
end

n=size(x,1);

fx=fft(x,[],1);

P2=abs(fx./L);
P1=P2(1:n/2+1,:);
P1(2:end-1,:)=2*P1(2:end-1,:);

f=Fs*((0:(n/2))/n)';
a=P1;

if ~isnan(frange)
    a=a(f>=frange(1) & f<=frange(2),:);
    f=f(f>=frange(1) & f<=frange(2));
end

if ~isnan(df)
    a=moving_smooth_tsm(a,f,df);
end

if ~isnan(figflag)
%     figure(figflag)
    if strcmp(plottype,'plot')==1
        plot(f,a)
    elseif strcmp(plottype,'semilog')==1
        semilog(f,a)
    else
        loglog(f,a)
    end
    xlabel('f / Hz')
    ylabel('Amplitude')
end

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

ax=zeros(size(a1));
for i=1:length(f)
    ax(i,:)=sqrt(mean((a1(f>(f(i).*10.^-df) & f<(f(i).*10.^+df),:)).^2,1));
    % log10(f)-log10(f1)=df
end

%a=interp1(f1,ax,f);

end
        