clear;
function [dydx]=der(y,x)
    n=size(y)
    dydx=diff(y,1)
    h=(max(x)-min(x))/size(x,1)
endfunction


hawaii=read('hawaii_melt_clap_impedance.dat',-1,3)
melt_hawaii=hawaii(:,1)
clap_hawaii=hawaii(:,2)
eta_hawaii=(hawaii(:,3))
n_hawaii=size(melt_hawaii,1)
m_hawaii=sqrt(n_hawaii)
dmeltdv_hawaii=zeros(m_hawaii,m_hawaii-1)
dmeltdclap_hawaii=zeros(m_hawaii,m_hawaii-1)
for ii=1:m_hawaii
   temp1_hawaii=diff(melt_hawaii((ii-1)*m_hawaii+1:ii*m_hawaii))./diff(eta_hawaii((ii-1)*m_hawaii+1:ii*m_hawaii))
   dmeltdv_hawaii(ii,:)=temp1_hawaii'
end
indx_hawaii=1:10:n_hawaii
for ii=1:m_hawaii
    temp1_hawaii=melt_hawaii(indx_hawaii+ii-1)
    temp2_hawaii=clap_hawaii(indx_hawaii+ii-1)
    temp3_hawaii=diff(temp1_hawaii)./diff(temp2_hawaii)
    dmeltdclap_hawaii(ii,:)=temp3_hawaii'
end
clapsd_hawaii=st_deviation(linspace(min(clap_hawaii),max(clap_hawaii),m_hawaii-1))

etasd_hawaii=st_deviation(linspace(min(eta_hawaii),max(eta_hawaii),m_hawaii-1))

for ii=1:m_hawaii
    err_hawaii(ii)=sqrt(clapsd_hawaii.^2.*(mean(dmeltdclap_hawaii(ii,:))).^2+etasd_hawaii.^2.*(mean(dmeltdv_hawaii(ii,:))).^2)
end

//dmeltdclap=der(melt,clap)
//dmeltdeta=der(melt,eta)
//



coral_sea=read('coral_sea_melt_clap_impedance.dat',-1,3)
melt_coral_sea=coral_sea(:,1)
clap_coral_sea=coral_sea(:,2)
eta_coral_sea=(coral_sea(:,3))
n_coral_sea=size(melt_coral_sea,1)
m_coral_sea=sqrt(n_coral_sea)
dmeltdv_coral_sea=zeros(m_coral_sea,m_coral_sea-1)
dmeltdclap_coral_sea=zeros(m_coral_sea,m_coral_sea-1)
for ii=1:m_coral_sea
   temp1_coral_sea=diff(melt_coral_sea((ii-1)*m_coral_sea+1:ii*m_coral_sea))./diff(eta_coral_sea((ii-1)*m_coral_sea+1:ii*m_coral_sea))
   dmeltdv_coral_sea(ii,:)=temp1_coral_sea'
end
indx_coral_sea=1:10:n_coral_sea
for ii=1:m_coral_sea
    temp1_coral_sea=melt_coral_sea(indx_coral_sea+ii-1)
    temp2_coral_sea=clap_coral_sea(indx_coral_sea+ii-1)
    temp3_coral_sea=diff(temp1_coral_sea)./diff(temp2_coral_sea)
    dmeltdclap_coral_sea(ii,:)=temp3_coral_sea'
end
clapsd_coral_sea=st_deviation(linspace(min(clap_coral_sea),max(clap_coral_sea),m_coral_sea-1))

etasd_coral_sea=st_deviation(linspace(min(eta_coral_sea),max(eta_coral_sea),m_coral_sea-1))

for ii=1:m_coral_sea
    err_coral_sea(ii)=sqrt(clapsd_coral_sea.^2.*(mean(dmeltdclap_coral_sea(ii,:))).^2+etasd_coral_sea.^2.*(mean(dmeltdv_coral_sea(ii,:))).^2)
end

//dmeltdclap=der(melt,clap)
//dmeltdeta=der(melt,eta)
//
xset("window",0)
clf
subplot(2,1,1)
plot2d(clap_hawaii,melt_hawaii,-11)
b=gca()
b.children(1).children.mark_size=2
b.children(1).children.mark_foreground=2
b.children(1).children.mark_background=0
plot2d(clap_coral_sea,melt_coral_sea,-9)
b=gca()
b.children(1).children.mark_size=2
b.children(1).children.mark_foreground=5
b.children(1).children.mark_background=0
subplot(2,1,2)
plot2d(eta_hawaii,melt_hawaii,-11)
b=gca()
b.children(1).children.mark_size=2
b.children(1).children.mark_foreground=2
b.children(1).children.mark_background=0
plot2d(eta_coral_sea,melt_coral_sea,-9)
b=gca()
b.children(1).children.mark_size=2
b.children(1).children.mark_foreground=5
b.children(1).children.mark_background=0

legend(['Hawaii' 'Coral Sea'],2,boxed=%f)
