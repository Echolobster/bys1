function RS2G_SimulatorParallelExperimentResults(varargin)
%% RS2G_SimulatorParallelExperimentResults
%
%   RS2G_SimulatorParallelExperimentResults()
%
%   Plots results multiple simulations using analyses developed to tease
%   apart models.
%
%
%%

%% Defaults
PlotOpts_default.colors = [0 0 1; 1 0 0];


%% Parse input

Parser = inputParser;

addParameter(Parser,'ModelSimulated','BLS')
addParameter(Parser,'Models',{'BLS','EKF','LNE','BLS_{mem}'})
addParameter(Parser,'resampleN',9)
addParameter(Parser,'replacement',false)
addParameter(Parser,'PlotOpts',PlotOpts_default)

parse(Parser,varargin{:})

ModelSimulated = Parser.Results.ModelSimulated;
Models = Parser.Results.Models;
resampleN = Parser.Results.resampleN;
replacement = Parser.Results.replacement;
PlotOpts = Parser.Results.PlotOpts;

%% Set up matrices
wmTrue = nan(100,1);
wpTrue = nan(100,2);
bTrue = nan(100,1);
lapseTrue = nan(100,1);
wm_dirftTrue = nan(100,1);

RMSE = nan(100,2);
BIAS = nan(100,2);
VAR = nan(100,2);

wmModels = nan(100,4);
wpModels = nan(100,4);
bModels = nan(100,4);
lapseModels = nan(100,4);
wm_driftModels = nan(100,4);
LLmodels = nan(100,4);

SimRMSE = nan(100,4,2);
SimVAR = nan(100,4,2);
SimBIAS = nan(100,4,2);

%% Get data
for simi = 1:100
    if exist([ModelSimulated '_All' num2str(simi) '.mat'],'file')
        load([ModelSimulated '_All' num2str(simi) '.mat'])
        
        modelTrue{simi} = model;
        wmTrue(simi,:) = wm;
        wpTrue(simi,:) = wp;
        bTrue(simi,:) = b;
        lapseTrue(simi,:) = lapse;
        wm_driftTrue(simi,:) = wm_drift;
        
        RMSE(simi,:) = rmseTrue(simi,:);
        BIAS(simi,:) = sqrt(biasTrue(simi,:));
        VAR(simi,:) = sqrt(variance(simi,:));
        
        wmModels(simi,:) = mean(FitResults.wm(simi,:,:),3);
        wpModels(simi,:) = mean(FitResults.wp(simi,:,:),3);
        bModels(simi,:) = mean(FitResults.b(simi,:,:),3);
        lapseModels(simi,:) = mean(FitResults.lapse(simi,:,:),3);
        wm_driftModels(simi,:) = mean(FitResults.wm_drift(simi,:,:),3);
        
        LLmodels(simi,:) = mean(FitResults.ll(simi,:,:),3);
        
        SimRMSE(simi,:,:) = FitResults.SimRMSE(simi,:,:);
        SimVAR(simi,:,:) = FitResults.SimVar(simi,:,:);
        SimBIAS(simi,:,:) = FitResults.SimBias(simi,:,:);
        
        TSs{simi} = TS;
        TPs{simi} = TP;
        MTP{simi} = mtp;
        STDtp{simi} = stdtp;
        TAs{simi} = ta;
    end
end

%% Resample
if resampleN
    Inds = find(~isnan(wmTrue(:,1)));
    sampleInds = Inds(randsample(length(Inds),resampleN));%,replacement));
    wmTrue = wmTrue(sampleInds,:);
    wpTrue = wpTrue(sampleInds,:);
    bTrue = bTrue(sampleInds,:);
    lapseTrue = lapseTrue(sampleInds,:);
    wm_driftTrue = wm_driftTrue(sampleInds,:);
    
    RMSE = RMSE(sampleInds,:);
    BIAS = BIAS(sampleInds,:);
    VAR = VAR(sampleInds,:);
    
    wmModels = wmModels(sampleInds,:);
    wpModels = wpModels(sampleInds,:);
    bModels = bModels(sampleInds,:);
    lapseModels = lapseModels(sampleInds,:);
    wm_driftModels = wm_driftModels(sampleInds,:);
    
    LLmodels = LLmodels(sampleInds,:);
    SimRMSE = SimRMSE(sampleInds,:,:);
    SimVAR = SimVAR(sampleInds,:,:);
    SimBIAS = SimBIAS(sampleInds,:,:);
    
    
    modelTrue = modelTrue(sampleInds);
    TSs = TSs(sampleInds);
    TPs = TPs(sampleInds);
    MTP = MTP(sampleInds);
    STDtp = STDtp(sampleInds);
    TAs = TAs(sampleInds);
else
    resampleN = 100;
end

modelInd = find(strcmp(Models,ModelSimulated));

%% Plot Results

%% Behavior/model fit
 %% TS v TP
 fH1 = figure('Name','t_s v t_p','Units','normalized','Position',[0.0201 0.0900 0.9632 0.8244]);
 sInd(1) = randsample(find(wmTrue == min(wmTrue)),1);
 sInd(2) = randsample(find(wmTrue == max(wmTrue)),1);
 for si = 1:2
     subplot(1,2,si)
     plotind = 1;
     ax = [400 1200 400 1200]; %[min(allts)-100 max(allts)+100 min(alltp)-100 max(alltp)+100];
     xticks = ax(1)+100:200:ax(2)-100;
     xticklabels = strread(num2str(xticks),'%s');
     yticks = ax(3)+100:200:ax(4)-100;
     yticklabels = strread(num2str(yticks),'%s');
     axis(ax);
     plotUnity;
     hold on
     
     for i = 1:N
         eh = errorbar(tss(:,1),MTP{sInd(si)}{i}(:,1),STDtp{sInd(si)}{i}(:,1),...
             'o','Color',PlotOpts.colors(i,:)+0.6*(~PlotOpts.colors(i,:)),'LineWidth',2);
         set(eh,'MarkerFaceColor',PlotOpts.colors(i,:))
         hold all
         plot(tsvec,TAs{sInd(si)}{i}(:,1,modelInd)+bModels(sInd(si),modelInd),...
             'Color',PlotOpts.colors(i,:)+0.6*(~PlotOpts.colors(i,:)),'LineWidth',2)
         plot(tss(:,1),MTP{sInd(si)}{i}(:,1),'o','Color',PlotOpts.colors(i,:),...
             'MarkerFaceColor',PlotOpts.colors(i,:),'MarkerSize',10)
         
         axis(ax)
         axis square
         xlabel('t_s(ms)')
         ylabel('t_p (ms)')
         mymakeaxis(gca,'xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels);
     end
 end

%% Parameters

figure('Name','Parameters','Position',[191 194 890 549])
subplot(2,3,1)
plot(wmTrue,wmModels(:,modelInd),'ko','MarkerFaceColor','k')
hold on
xlabel('w_m simulated')
ylabel('w_m fit')
axis([0.01 0.3 0.01 0.3])
plotUnity;
axis square
mymakeaxis(gca)

subplot(2,3,2)
plot(wpTrue(:,1),wpModels(:,modelInd),'ko','MarkerFaceColor','k')
hold on
xlabel('w_p simulated')
ylabel('w_p fit')
axis([0.01 0.2 0.01 0.2])
plotUnity;
axis square
mymakeaxis(gca)

subplot(2,3,3)
plot(bTrue(:,1),bModels(:,modelInd),'ko','MarkerFaceColor','k')
hold on
xlabel('b simulated')
ylabel('b fit')
axis([0 50 0 50])
plotUnity;
axis square
mymakeaxis(gca)

subplot(2,3,4)
plot(lapseTrue(:,1),lapseModels(:,modelInd),'ko','MarkerFaceColor','k')
hold on
xlabel('\lambda simulated')
ylabel('\lambda fit')
axis([0 1 0 1])
plotUnity;
axis square
mymakeaxis(gca)

subplot(2,3,5)
plot(wm_driftTrue(:,1),wm_driftModels(:,modelInd),'ko','MarkerFaceColor','k')
hold on
xlabel('w_{mem} simulated')
ylabel('w_{mem} fit')
axis([0.01 0.4 0.01 0.4])
plotUnity;
axis square
mymakeaxis(gca)

%% LL models

figure('Name','Log Likelihood Models')
plot((LLmodels - repmat(LLmodels(:,1),[1,size(LLmodels,2)]))','ko',...
    'MarkerFaceColor','k');
xticklabels = {'BLS','EKF','LNE','BLS_{mem}'};
xlabel('Model')
ylabel('Relative LL (compared to BLS)')
hold on
plotHorizontal(0);
mymakeaxis(gca,'xticks',[1 2 3 4],'xticklabels',xticklabels)

%% Normalized RMSE
del = 0.1;
fh = figure('Name','normalized RMSE','Position', [150 35 1036 763]);
for modeli = 1:length(Models)
    subplot(2,2,modeli)
    plot([1 2],[1 1],'k')
    hold on
    for i = 1:resampleN
        h(i) = plot([1 2],RMSE(i,:)./repmat(SimRMSE(i,modeli,2),1,2),...
            'k-o');
        set(h(i),'MarkerFaceColor',h(i).Color)
    end
    errorbar([1-del 2+del],nanmean(RMSE./repmat(SimRMSE(:,modeli,2),1,2),1),...
        nanstd(RMSE./repmat(SimRMSE(:,modeli,2),1,2),[],1)/sqrt(sum(~isnan(wmTrue(:,1)))),...
        'ko','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','k')
    axis square
    xlabel('Condition')
    ylabel('RMSE_x/expectedRMSE_2')
    T{modeli} = ['Normalized to ' Models{modeli} ' model expectation'];
    text(1,1.002,'Optimal performance')
    
    ax(modeli,:) = axis;
end

for modeli = 1:length(Models)
    subplot(2,2,modeli)
    axis([min(ax(:,1)) max(ax(:,2)) min(ax(:,3)) max(ax(:,4))])
    ah = gca;
    xticks = [1 2];
    xticklabels = {'RSG','RSSG'};
    yticks = ah.YTick(1:2:end);
    yticklabels = ah.YTickLabel(1:2:end);
    mymakeaxis(gca,'xytitle',T{modeli},'xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
end

%% deltaRMSE
figure('Name','\DeltaRMSE','Position', [150 35 1036 763])
edgesM = linspace(-12,12,21);
for modeli = 1:4
    [~,deltaP(modeli),~,deltaSTATS(modeli)] = ...
        ttest((RMSE(:,2)-RMSE(:,1))-(SimRMSE(:,modeli,2)-SimRMSE(:,modeli,1)));
    
    subplot(2,2,modeli)
    [n, edges, ...
     dataH, patchH, zAxH, nAxH, ztextH, ntextH, zAxTickH, nAxTickH, unityH] = ...
     diagonalMarginal(RMSE(:,2)-RMSE(:,1),SimRMSE(:,modeli,2)-SimRMSE(:,modeli,1),...
                'binN',9,'T',[0.1 0.1],'nscale',8,'edges',edgesM,'offsetFactor',1,...
                'axVals',[-20 10 -20 10]);
end

for modeli = 1:length(Models)
    subplot(2,2,modeli)
    text(9,9,['t(' num2str(deltaSTATS(modeli).df) ') = ' num2str(deltaSTATS(modeli).tstat)])
    text(8,8,['p-value = ' num2str(deltaP(modeli))])
    xlabel('\Delta RMSE_{Observed} (ms)')
    ylabel(['\Delta RMSE_{' Models{modeli} '} (ms)'])
    ah = gca;
    xticks = -20:10:0;
    xticklabels = {'-20','-10','0'};
    yticks = -20:10:0;
    yticklabels = {'-20','-10','0'};
    mymakeaxis(ah,'xytitle',Models{modeli},...
        'xticks',xticks,'xticklabels',xticklabels,'yticks',yticks,'yticklabels',yticklabels)
end

%% Residual BIAS
figure('Name','Residual BIAS','Position',[155 296 1212 372])
for modeli = 1:length(Models)
    subplot(2,2,modeli)
    residualBias(:,:,modeli) = BIAS-permute(SimBIAS(:,modeli,:),[1 3 2]);
    for j = 1:2
        [~, pvalBias(:,j,modeli), ~, STATSbias(:,j,modeli)]= ttest(residualBias(:,j,modeli));
    end
    edges = linspace(-10,10,8);
    x = edges + (edges(2) - edges(1))/2;
    ncount = histcounts(residualBias(:,1,modeli),edges);
    barProperties.FaceColor = PlotOpts.colors(1,:);
    barProperties.EdgeColor = 'none';
    barProperties.ShowBaseLine = 'off';
    barProperties.BarWidth = 0.8;
    barProperties.FaceAlpha = 0.3;
    h = mybargraph(x(1:end-1),ncount,'barProperties',barProperties);
    hold on
    
    text(-8,2,['p-val = ' num2str(pvalBias(:,1,modeli))],'Color',PlotOpts.colors(1,:));
    
    ncount = histcounts(residualBias(:,2,modeli),edges);
    barProperties.FaceColor = PlotOpts.colors(2,:);
    barProperties.EdgeColor = 'none';
    barProperties.ShowBaseLine = 'off';
    barProperties.BarWidth = 0.8;
    barProperties.FaceAlpha = 0.3;
    mybargraph(x(1:end-1),ncount,'barProperties',barProperties);
    
    text(-8,2.5,['p-val = ' num2str(pvalBias(:,2,modeli))],'Color',PlotOpts.colors(2,:));
    
    %axis([-10 10 0 6])
    
    ylabel('Count')
    xlabel('Residual BIAS (ms)')
    mymakeaxis(gca,'xytitle',Models{modeli},...
        'xticks',[-10 0 10],'xticklabels',{'-10','0','10'});
end

%% Residual VAR
figure('Name','Residual VAR','Position',[155 296 1212 372])
for modeli = 1:length(Models)
    subplot(2,2,modeli)
    residualVar(:,:,modeli) = VAR-permute(SimVAR(:,modeli,:),[1 3 2]);
    for j = 1:2
        [~, pvalVar(:,j,modeli), ~, STATSvar(:,j,modeli)]= ttest(residualVar(:,j,modeli));
    end
    edges = linspace(-10,10,8);
    x = edges + (edges(2) - edges(1))/2;
    ncount = histcounts(residualVar(:,1,modeli),edges);
    barProperties.FaceColor = PlotOpts.colors(1,:);
    barProperties.EdgeColor = 'none';
    barProperties.ShowBaseLine = 'off';
    barProperties.BarWidth = 0.8;
    barProperties.FaceAlpha = 0.3;
    h = mybargraph(x(1:end-1),ncount,'barProperties',barProperties);
    hold on
    
    text(-8,4,['p-val = ' num2str(pvalVar(:,1,modeli))],'Color',PlotOpts.colors(1,:));
    
    ncount = histcounts(residualVar(:,2,modeli),edges);
    barProperties.FaceColor = PlotOpts.colors(2,:);
    barProperties.EdgeColor = 'none';
    barProperties.ShowBaseLine = 'off';
    barProperties.BarWidth = 0.8;
    barProperties.FaceAlpha = 0.3;
    mybargraph(x(1:end-1),ncount,'barProperties',barProperties);
    
    text(-8,4.5,['p-val = ' num2str(pvalVar(:,2,modeli))],'Color',PlotOpts.colors(2,:));
    
    %axis([-10 10 0 6])
    
    ylabel('Count')
    xlabel('Residual VAR (ms)')
    mymakeaxis(gca,'xytitle',Models{modeli},...
        'xticks',[-10 0 10],'xticklabels',{'-10','0','10'});
end