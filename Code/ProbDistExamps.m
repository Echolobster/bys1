%% ProbDistExamps
%
%   Creates examples likelihoods, prior and posterior distributions for the
%   RSNG experiments
%
%%

% Variables
N = 2;              % Number of measurements
tsmin = 600;
tsmax = 1000;
examplets = [800];
exampletms = [700 900];
exampleN = [1 2];
dts = 0.1;
tsvec = 300:dts:1400;
tms = 200:1600;

wm = 0.15;  % scalar variance of measurements; conditionally independent for now
wp = 0.1;


% Likelihood funtion
Likelihoodf = @(w,N,tm,ts)( (1./(sqrt(2*pi)*w*ts)).^N .* exp( -(tm-ts)'*(tm-ts) ./ (2*w.^2.*ts.^2) ) );

% Eqivalent likelihood for N measurements
eLF = @(w,N,tm,tmeff,ts)( (1./(sqrt(2*pi)*w*ts)).^N .* exp( (N*tmeff^2)/(2*w^2*ts.^2) ) .* exp( -sum(tm.^2)./(2*w^2*ts.^2) ) .* exp( -(tmeff-ts).^2/(2*(w*ts/sqrt(N))^2) ) );

% Generate noisy measurements
for i = 1:length(examplets)
%     noise = wm*(examplets(i)*ones(1,N)).*randn(1,N);
%     tm(i,:) = examplets(i)*ones(1,N) + noise;
    tm(i,:) = exampletms;
end

% Generate distributions
%Prior
p0 = zeros(size(tsvec));
p0(find(tsvec == tsmin):find(tsvec == tsmax)) = ones(size(tsmin:dts:tsmax))./(tsmax - tsmin)*dts;
inds = 0;
for nn = exampleN
    inds = inds+1;
    for i = 1:length(examplets)        
        % Likelihoods
        for n = 1:nn
            l(:,n,i,nn) = arrayfun(@(tsi)(Likelihoodf(wm,1,tm(i,n),tsi)),tsvec);
        end
        
        L(:,i,inds) = arrayfun(@(tsi)(Likelihoodf(wm,nn,tm(i,1:nn)',tsi)),tsvec)/sum((arrayfun(@(tsi)(Likelihoodf(wm,nn,tm(i,1:nn)',tsi)),tsvec)));
        eL(:,i,inds) = arrayfun(@(tsi)(eLF(wm,nn,tm(i,1:nn),mean(tm(i,1:nn)),tsi)),tsvec)/sum((arrayfun(@(tsi)(eLF(wm,nn,tm(i,1:nn),mean(tm(i,1:nn)),tsi)),tsvec)));
        
        % Posterior
        posterior(:,i,inds) = (arrayfun(@(tsi)(Likelihoodf(wm,nn,tm(i,1:nn)',tsi)),tsvec).*p0)/sum((arrayfun(@(tsi)(Likelihoodf(wm,nn,tm(i,1:nn)',tsi)),tsvec)).*p0);
        
        te(inds,i) = sum(tsvec'.*posterior(:,i,inds));
    end
end

% Building of likelihood function
P(:,1) = arrayfun(@(tmi)(Likelihoodf(wm,1,tmi,tsmin)),tms);
P(:,2) = arrayfun(@(tmi)(Likelihoodf(wm,1,tmi,(tsmin+tsmax)/2)),tms);
P(:,3) = arrayfun(@(tmi)(Likelihoodf(wm,1,tmi,tsmax)),tms);

lf = arrayfun(@(tsi)(Likelihoodf(wm,1,(tsmin+tsmax)/2,tsi)),tms);

%% Plot
for nn= 1:length(exampleN)
    h(nn) = figure;
    % Prior
    subplot(3,1,2)
    plot(tsvec,p0,'Color',[0.32 0.32 0.32],'LineWidth',2);
    xlabel('t_s (ms)')
    ylabel('Probability')
    axis tight
    ax2 = axis;
    ax1 = axis;
    ax3 = axis;
    ax3(4) = 0;
    for i = 1:length(examplets)
        if i == 1
            m = 0;
        else
            m = 1;
        end
        % Likelihood
        subplot(3,1,3)
        ltemp = squeeze(l(:,1:exampleN(nn),i,exampleN(nn)));
        plot(tsvec,ltemp./repmat(sum(ltemp,1),size(ltemp,1),1),'Color',[0.64 0.64 0.64]-m*[0.64 0.64 0.64])
        hold on
        plot(tsvec,L(:,i,nn),'Color',[0.64 0.64 0.64]-m*[0.64 0.64 0.64],'LineWidth',2)
        axis tight
        temp = axis;
        ax3(temp > ax3) = temp(temp > ax3);
        xlabel('t_s (ms)')
        ylabel('Relative probability')
        
        % Posterior
        subplot(3,1,1)
        plot(tsvec,posterior(:,i,nn),'Color',[0.64 0.64 0.64]-m*[0.64 0.64 0.64],'LineWidth',2);
        hold on
        plot(ceil(te(nn,i)),posterior(tsvec == ceil(te(nn,i)),i,nn),'k.','MarkerSize',20)
        axis tight
        temp = axis;
        ax1(temp > ax1) = temp(temp > ax1);
        xlabel('t_s (ms)')
        ylabel('Probability')
    end
%     axmin = min([ax1(3) ax2(3) ax3(3)]);
%     axmax = max([ax1(4) ax2(4) ax3(4)]);
%     ax1(3) = axmin;
%     ax2(3) = axmin;
%     ax3(3) = axmin;
%     ax1(4) = axmax;
%     ax2(4) = axmax;
%     ax3(4) = axmax;
end

for nn = 1:length(exampleN)
    figure(h(nn))
    for i = 1:length(examplets)
        if i == 1
            m = 0;
        else
            m = 1;
        end
        subplot(3,1,3)
        ltemp = squeeze(l(:,1:exampleN(nn),i,exampleN(nn)));
        q = find(ltemp == repmat(max(ltemp),size(ltemp,1),1));
        for j = 1:length(q)
            q(j) = q(j) - (j-1)*length(tsvec);
            plot(tsvec(q(j))*[1 1],[0 ax3(4)],'--','Color',[0.64 0.64 0.64]-m*[0.64 0.64 0.64])
        end
        
        plot(((-1+sqrt(1+4*wm^2))/(2*wm^2))*mean(tm(i,1:exampleN(nn)))*[1 1],[0 ax3(4)],'-.','Color',[0.64 0.64 0.64]-m*[0.64 0.64 0.64])
        plot(tsvec(L(:,i,nn) == max(L(:,i,nn)))*[1 1],[0 ax3(4)],'Color',[0.64 0.64 0.64]-m*[0.64 0.64 0.64])
    end
end

ax2(4) = 2*ax2(4);
for nn = 1:length(exampleN)
    figure(h(nn))
    subplot(3,1,1)
    axis(ax1)
    xytitle = 'Posterior';
    mymakeaxis(gca,'xytitle',xytitle,'xticks',400:200:1400,'xticklabels',strread(num2str(400:200:1400),'%s'),'yticks',[0:2:4]*10^-4,'yticklabels',strread(num2str([0:2:4]*10^-4),'%s'))
    subplot(3,1,2)
    axis(ax2)
    xytitle = 'Posterior';
    mymakeaxis(gca,'xytitle',xytitle,'xticks',400:200:1400,'xticklabels',strread(num2str(400:200:1400),'%s'),'yticks',[0:2:4]*10^-4,'yticklabels',strread(num2str([0:2:4]*10^-4),'%s'))
    subplot(3,1,3)
    axis(ax3)
    xytitle = 'Posterior';
    mymakeaxis(gca,'xytitle',xytitle,'xticks',400:200:1400,'xticklabels',strread(num2str(400:200:1400),'%s'),'yticks',[0:2:4]*10^-4,'yticklabels',strread(num2str([0:2:4]*10^-4),'%s'))
end
    

figure
subplot(2,1,1)
for i = 1:3
    plot(tms,P(:,i),'Color',[0.66 0.66 0.66] - (i-1)*[0.33 0.33 0.33],'LineWidth',2)
    hold on
end
ax = axis;
plot((tsmin+tsmax)/2*[1 1],[0 ax(4)],'r','LineWidth',2)
xlabel('t_m (ms)')
ylabel('P(t_m|t_s)')
legend(['t_s = ' num2str(tsmin)],['t_s = ' num2str((tsmin+tsmax)/2)],['t_s = ' num2str(tsmax)],['t_m = ' num2str((tsmin+tsmax)/2)])
mymakeaxis(gca,'xticks',200:400:1400,'xticklabels',strread(num2str(200:400:1400),'%s'),'yticks',(0:2:4)*10^-3,'yticklabels',strread(num2str((0:2:4)*10^-3),'%s'))
subplot(2,1,2)
plot(tms,lf,'r','LineWidth',2)
hold on
plot(tsmin,lf(tms == tsmin),'.','MarkerSize',20,'Color',[0.66 0.66 0.66])
plot((tsmin+tsmax)/2,lf(tms == (tsmin+tsmax)/2),'.','MarkerSize',20,'Color',[0.33 0.33 0.33])
plot(tsmax,lf(tms == tsmax),'.','MarkerSize',20,'Color',[0 0 0])
xlabel('t_s (ms)')
ylabel('Likelihood')
mymakeaxis(gca,'xticks',200:400:1400,'xticklabels',strread(num2str(200:400:1400),'%s'),'yticks',(0:1:3)*10^-3,'yticklabels',strread(num2str((0:1:3)*10^-3),'%s'))

%% BLS model
v  = 700:50:900;
method_opts.dx = 10;
method_opts.type = 'quad';
tms = 400:1200;
[TM1, TM2] = meshgrid(tms);
M = [TM1(:) TM2(:)];
fBLS2 = ScalarBayesEstimators(M,wm,tsmin,tsmax,'method',method_opts);
figure('Name','fBLS2')
imagesc(tms,tms,reshape(fBLS2,length(tms),length(tms)))
axis xy
colormap gray
axis square
a = [595 1010];%caxis;
hold on
colormap gray
[C, H] = contour(TM1,TM2,reshape(fBLS2,length(tms),length(tms)),v,'Color',[1 215/255 0],'LineWidth',1.5);
caxis(a);
xlabel('t_{m_1}')
ylabel('t_{m_2}')
mymakeaxis(gca,'xticks',400:200:1200,'xticklabels',strread(num2str(400:200:1200),'%s'),'yticks',400:200:1200,'yticklabels',strread(num2str(400:200:1200),'%s'))
ax = gca;
h = colorbar;
h.TickDirection = 'out';
h.Box = 'off';
h.Position(1) = h.Position(1)+0.1;
h.Position(2) = ax.Position(2) + 0.25;
h.Position(4) = ax.Position(4) - 0.25;

%% BLS update
fBLS1 = ScalarBayesEstimators(tms(:),wm,tsmin,tsmax,'method',method_opts);
figure('Name','fBLS update')
subplot(1,2,1)
imagesc(tms,tms,reshape(fBLS2,length(tms),length(tms))-repmat(fBLS1,[1,size(TM2,2)])')
axis xy
colormap gray
axis square
a2 = caxis;
hold on
colormap gray
[C, H] = contour(TM1,TM2,reshape(fBLS2,length(tms),length(tms))-repmat(fBLS1,[1,size(TM2,2)])',-150:50:150,'Color',[1 215/255 0],'LineWidth',1.5);
caxis(a2);
xlabel('t_{m_1}')
ylabel('t_{m_2}')
mymakeaxis(gca,'xticks',400:200:1200,'xticklabels',strread(num2str(400:200:1200),'%s'),'yticks',400:200:1200,'yticklabels',strread(num2str(400:200:1200),'%s'))
ax = gca;
h = colorbar;
h.TickDirection = 'out';
h.Box = 'off';
h.Position(1) = h.Position(1)+0.1;
h.Position(2) = ax.Position(2) + 0.25;
h.Position(4) = ax.Position(4) - 0.25;

subplot(1,2,2)
imagesc(tms,tms,TM2-repmat(fBLS1,[1,size(TM2,2)])')
axis xy
colormap gray
axis square
a2 = caxis;
hold on
colormap gray
[C, H] = contour(TM1,TM2,TM2-repmat(fBLS1,[1,size(TM2,2)])',-300:100:300,'Color',[1 215/255 0],'LineWidth',1.5);
caxis(a2);
xlabel('t_{m_1}')
ylabel('t_{m_2}')
mymakeaxis(gca,'xticks',400:200:1200,'xticklabels',strread(num2str(400:200:1200),'%s'),'yticks',400:200:1200,'yticklabels',strread(num2str(400:200:1200),'%s'))
ax = gca;
h = colorbar;
h.TickDirection = 'out';
h.Box = 'off';
h.Position(1) = h.Position(1)+0.1;
h.Position(2) = ax.Position(2) + 0.25;
h.Position(4) = ax.Position(4) - 0.25;

%% EKF model
estimator.type = 'EKF';
fEKF2 = ScalarBayesEstimators(M,wm,tsmin,tsmax,'method',method_opts,'estimator',estimator);
figure('Name','fEKF2')
imagesc(tms,tms,reshape(fEKF2,length(tms),length(tms)))
axis xy
colormap gray
axis square
% a = caxis;
hold on
colormap gray
[C, H] = contour(TM1,TM2,reshape(fEKF2,length(tms),length(tms)),v,'Color',[1 215/255 0],'LineWidth',1.5);
caxis(a);
xlabel('t_{m_1}')
ylabel('t_{m_2}')
mymakeaxis(gca,'xticks',400:200:1200,'xticklabels',strread(num2str(400:200:1200),'%s'),'yticks',400:200:1200,'yticklabels',strread(num2str(400:200:1200),'%s'))
ax = gca;
h = colorbar;
h.TickDirection = 'out';
h.Box = 'off';
h.Position(1) = h.Position(1)+0.1;
h.Position(2) = ax.Position(2) + 0.25;
h.Position(4) = ax.Position(4) - 0.25;

%% EKF update
fBLS1 = ScalarBayesEstimators(tms(:),wm,tsmin,tsmax,'method',method_opts);
figure('Name','fBLS update')
subplot(1,2,1)
imagesc(tms,tms,reshape(fEKF2,length(tms),length(tms))-repmat(fBLS1,[1,size(TM2,2)])')
axis xy
colormap gray
axis square
a2 = caxis;
hold on
colormap gray
[C, H] = contour(TM1,TM2,reshape(fEKF2,length(tms),length(tms))-repmat(fBLS1,[1,size(TM2,2)])',-150:50:150,'Color',[1 215/255 0],'LineWidth',1.5);
caxis(a2);
xlabel('t_{m_1}')
ylabel('t_{m_2}')
mymakeaxis(gca,'xticks',400:200:1200,'xticklabels',strread(num2str(400:200:1200),'%s'),'yticks',400:200:1200,'yticklabels',strread(num2str(400:200:1200),'%s'))
ax = gca;
h = colorbar;
h.TickDirection = 'out';
h.Box = 'off';
h.Position(1) = h.Position(1)+0.1;
h.Position(2) = ax.Position(2) + 0.25;
h.Position(4) = ax.Position(4) - 0.25;

subplot(1,2,2)
imagesc(tms,tms,TM2-repmat(fBLS1,[1,size(TM2,2)])')
axis xy
colormap gray
axis square
a2 = caxis;
hold on
colormap gray
[C, H] = contour(TM1,TM2,TM2-repmat(fBLS1,[1,size(TM2,2)])',-300:100:300,'Color',[1 215/255 0],'LineWidth',1.5);
caxis(a2);
xlabel('t_{m_1}')
ylabel('t_{m_2}')
mymakeaxis(gca,'xticks',400:200:1200,'xticklabels',strread(num2str(400:200:1200),'%s'),'yticks',400:200:1200,'yticklabels',strread(num2str(400:200:1200),'%s'))
ax = gca;
h = colorbar;
h.TickDirection = 'out';
h.Box = 'off';
h.Position(1) = h.Position(1)+0.1;
h.Position(2) = ax.Position(2) + 0.25;
h.Position(4) = ax.Position(4) - 0.25;

%% BLS_wm1wm2
estimator.type = 'BLS_wm1wm2';
estimator.wm_drift = 1.5*wm;
E = ScalarBayesEstimators(M,wm,tsmin,tsmax,'method',method_opts,'estimator',estimator);
figure('Name','fBLS_wm1wm2')
imagesc(tms,tms,reshape(E,length(tms),length(tms)))
axis xy
colormap gray
axis square
% a = caxis;
hold on
colormap gray
[C, H] = contour(TM1,TM2,reshape(E,length(tms),length(tms)),v,'Color',[1 215/255 0],'LineWidth',1.5);
caxis(a);
xlabel('t_{m_1}')
ylabel('t_{m_2}')
mymakeaxis(gca,'xticks',400:200:1200,'xticklabels',strread(num2str(400:200:1200),'%s'),'yticks',400:200:1200,'yticklabels',strread(num2str(400:200:1200),'%s'))
ax = gca;
h = colorbar;
h.TickDirection = 'out';
h.Box = 'off';
h.Position(1) = h.Position(1)+0.1;
h.Position(2) = ax.Position(2) + 0.25;
h.Position(4) = ax.Position(4) - 0.25;

%% Linear-nonlienar estimator (LNE)
estimator.type = 'weightedMean';
estimator.weights = ones(1,n)/n;
E = ScalarBayesEstimators(M,wm,tsmin,tsmax,'method',method_opts,'estimator',estimator);
figure('Name','Averaging')
imagesc(tms,tms,reshape(E,length(tms),length(tms)))
axis xy
colormap gray
axis square
hold on
colormap gray
[C, H] = contour(TM1,TM2,reshape(E,length(tms),length(tms)),v,'Color',[1 215/255 0],'LineWidth',1.5);
% clabel(C,H)
caxis(a);
xlabel('t_{m_1}')
ylabel('t_{m_2}')
mymakeaxis(gca,'xticks',400:200:1200,'xticklabels',strread(num2str(400:200:1200),'%s'),'yticks',400:200:1200,'yticklabels',strread(num2str(400:200:1200),'%s'))
ax = gca;
h = colorbar;
h.TickDirection = 'out';
h.Box = 'off';
h.Position(1) = h.Position(1)+0.1;
h.Position(2) = ax.Position(2) + 0.25;
h.Position(4) = ax.Position(4) - 0.25;

%% MLE
estimator.type = 'MLE';
estimator.weights = ones(1,n)/n;
method_analytical.type = 'analytical';
MLE = ScalarBayesEstimators(M,wm,tsmin,tsmax,'method',method_analytical,'estimator',estimator);
figure('Name','MLE')
imagesc(tms,tms,reshape(MLE,length(tms),length(tms)))
axis xy
colormap gray
axis square
hold on
colormap gray
[C, H] = contour(TM1,TM2,reshape(MLE,length(tms),length(tms)),700:100:900,'Color',[1 215/255 0],'LineWidth',1.5);
%caxis(a);
xlabel('t_{m_1}')
ylabel('t_{m_2}')
mymakeaxis(gca,'xticks',400:200:1200,'xticklabels',strread(num2str(400:200:1200),'%s'),'yticks',400:200:1200,'yticklabels',strread(num2str(400:200:1200),'%s'))
ax = gca;
h = colorbar;
h.TickDirection = 'out';
h.Box = 'off';
h.Position(1) = h.Position(1)+0.1;
h.Position(2) = ax.Position(2) + 0.25;
h.Position(4) = ax.Position(4) - 0.25;


%% Average of Estimates
estimator.type = 'aveEstimates';
estimator.weights = ones(1,n)/n;
E = ScalarBayesEstimators(M,wm,tsmin,tsmax,'method',method_opts,'estimator',estimator);
figure('Name','Average of estimates')
imagesc(tms,tms,reshape(E,length(tms),length(tms)))
axis xy
colormap gray
axis square
hold on
colormap gray
[C, H] = contour(TM1,TM2,reshape(E,length(tms),length(tms)),v,'Color',[1 215/255 0],'LineWidth',1.5);
caxis(a);
xlabel('t_{m_1}')
ylabel('t_{m_2}')
mymakeaxis(gca,'xticks',400:200:1200,'xticklabels',strread(num2str(400:200:1200),'%s'),'yticks',400:200:1200,'yticklabels',strread(num2str(400:200:1200),'%s'))
ax = gca;
h = colorbar;
h.TickDirection = 'out';
h.Box = 'off';
h.Position(1) = h.Position(1)+0.1;
h.Position(2) = ax.Position(2) + 0.25;
h.Position(4) = ax.Position(4) - 0.25;

%% Single Measurement
fBLS1 = ScalarBayesEstimators(tms(:),wm,tsmin,tsmax,'method',method_opts);
figure('Name','fBLS1')
plot(tms,fBLS1,'k','LineWidth',2)
hold on
axis square
axis equal
plotUnity;
xlabel('t_{m}')
ylabel('t_a')
mymakeaxis(gca,'xticks',400:200:1200,'xticklabels',strread(num2str(400:200:1200),'%s'),'yticks',600:200:1000,'yticklabels',strread(num2str(600:200:1000),'%s'))

%% p(tm1,tm2,|ts,tp)
tms2 = 400:5:1200;
[TM12, TM22] = meshgrid(tms2);
M2 = [TM12(:) TM22(:)];
ptm1tm2_take_tstp = ScalarMeasuremnetProbabilityGivenInputOutput(800,800,M2,wm,wp,tsmin,tsmax);
figure('Name','Probability of t_m1 and t_m2, given ts and tp','Position',[680 270 1026 828])
imagesc(tms2,tms2,reshape(ptm1tm2_take_tstp,length(tms2),length(tms2)))
axis xy
axis square
hold on
colormap gray
cax = caxis;
[C, H] = contour(TM1,TM2,reshape(fBLS2,length(tms),length(tms)),700:100:900,'Color',[1 215/255 0],'LineWidth',1.5);
xlabel('\fontsize{16}{0}\selectfont $t_{m_1}$')
ylabel('\fontsize{16}{0}\selectfont $t_{m_2}$')
mymakeaxis(gca,'xticks',400:200:1200,...
    'xticklabels',strread(num2str(400:200:1200),'%s'),...
    'yticks',400:200:1200,'yticklabels',strread(num2str(400:200:1200),'%s'),...
    'interpreter','latex','xytitle','\fontsize{16}{0}\selectfont $p(t_{m_1},t_{m_2}|t_s = 800,t_p=800)$')
ax = gca;
caxis(cax);
h = colorbar;
h.TickDirection = 'out';
h.Box = 'off';
h.Position(1) = h.Position(1)+0.1;
h.Position(2) = ax.Position(2) + 0.25;
h.Position(4) = ax.Position(4) - 0.25;

%% p(tm1,tm2,|ts)
tms2 = 400:5:1200;
[TM12, TM22] = meshgrid(tms2);
M2 = [TM12(:) TM22(:)];
ptm1tm2_take_tstp = ScalarMeasuremnetProbabilityGivenInputOutput(800,NaN,M2,wm,wp,tsmin,tsmax);
figure('Name','Probability of t_m1 and t_m2, given ts and tp','Position',[680 270 1026 828])
imagesc(tms2,tms2,reshape(ptm1tm2_take_tstp,length(tms2),length(tms2)))
axis xy
axis square
hold on
colormap gray
cax = caxis;
[C, H] = contour(TM1,TM2,reshape(fBLS2,length(tms),length(tms)),700:100:900,'Color',[1 215/255 0],'LineWidth',1.5);
xlabel('\fontsize{16}{0}\selectfont $t_{m_1}$')
ylabel('\fontsize{16}{0}\selectfont $t_{m_2}$')
mymakeaxis(gca,'xticks',400:200:1200,...
    'xticklabels',strread(num2str(400:200:1200),'%s'),...
    'yticks',400:200:1200,'yticklabels',strread(num2str(400:200:1200),'%s'),...
    'interpreter','latex','xytitle','\fontsize{16}{0}\selectfont $p(t_{m_1},t_{m_2}|t_s = 800)$')
ax = gca;
caxis(cax);
h = colorbar;
h.TickDirection = 'out';
h.Box = 'off';
h.Position(1) = h.Position(1)+0.1;
h.Position(2) = ax.Position(2) + 0.25;
h.Position(4) = ax.Position(4) - 0.25;

%% p(tm1,tm2,|ts) for LNE model
estimator.type = 'weightedMean';
estimator.weights = ones(1,n)/n;
E = ScalarBayesEstimators(M,wm,tsmin,tsmax,'method',method_opts,'estimator',estimator);
tms2 = 400:5:1200;
[TM12, TM22] = meshgrid(tms2);
M2 = [TM12(:) TM22(:)];
ptm1tm2_take_tstp = ScalarMeasuremnetProbabilityGivenInputOutput(800,NaN,M2,wm,wp,tsmin,tsmax);
figure('Name','Probability of t_m1 and t_m2, given ts and tp, LNE','Position',[680 270 1026 828])
imagesc(tms2,tms2,reshape(ptm1tm2_take_tstp,length(tms2),length(tms2)))
axis xy
axis square
hold on
colormap gray
cax = caxis;
[C, H] = contour(TM1,TM2,reshape(E,length(tms),length(tms)),700:100:900,'Color',[1 215/255 0],'LineWidth',1.5);
xlabel('\fontsize{16}{0}\selectfont $t_{m_1}$')
ylabel('\fontsize{16}{0}\selectfont $t_{m_2}$')
mymakeaxis(gca,'xticks',400:200:1200,...
    'xticklabels',strread(num2str(400:200:1200),'%s'),...
    'yticks',400:200:1200,'yticklabels',strread(num2str(400:200:1200),'%s'),...
    'interpreter','latex','xytitle','\fontsize{16}{0}\selectfont $p(t_{m_1},t_{m_2}|t_s = 800)$')
ax = gca;
caxis(cax);
h = colorbar;
h.TickDirection = 'out';
h.Box = 'off';
h.Position(1) = h.Position(1)+0.1;
h.Position(2) = ax.Position(2) + 0.25;
h.Position(4) = ax.Position(4) - 0.25;