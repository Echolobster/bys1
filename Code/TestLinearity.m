function [pvalGreater, pvalLesser] = TestLinearity(slist,varargin)
%% TestLinearity
%
%   h = TestLinearity(slist)
%
%   Tests the hypothesis that reproduction data is linearly related to
%   sample interval.
%
%%

%% Defaults
TAexpectation_default.methodopts.dx = 0.01;
TAexpectation_default.dt = 10;
PlotOpts_default.colors = [0 0 1; 1 0 0];
TheoreticalRMSE_default.wmvec = NaN;
TheoreticalRMSE_default.type = 'EachSubject';

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'slist')     % List of subjects to analyze
addParameter(Parser,'N',2)      % Maximum number of sets
addParameter(Parser,'tss',600:100:1000)     % sample times for experiment
addParameter(Parser,'simulationN',10000)    % Number of trials per simulation
addParameter(Parser,'CommonFileName','_BLSbiasedFitResults20150913')
addParameter(Parser,'TAexpectation',TAexpectation_default)  % For controlling the calculation of the expected value of aim times under a model
addParameter(Parser,'TheoreticalRMSE',TheoreticalRMSE_default)
addParameter(Parser,'Plot','Yes')
addParameter(Parser,'PlotOpts',PlotOpts_default)
addParameter(Parser,'ExampleSubject',5)
addParameter(Parser,'alpha',0.05)           % Significance level

parse(Parser,slist,varargin{:})

slist = Parser.Results.slist;
N = Parser.Results.N;
tss = Parser.Results.tss;
simulationN = Parser.Results.simulationN;
CommonFileName = Parser.Results.CommonFileName;
TAexpectation = Parser.Results.TAexpectation;
TheoreticalRMSE = Parser.Results.TheoreticalRMSE;
Plot = Parser.Results.Plot;
PlotOpts = Parser.Results.PlotOpts;
ExampleSubject = Parser.Results.ExampleSubject;
alpha = Parser.Results.alpha;

%% Load data, regress, and find residuals, test if significantly different
for i = 1:length(slist)
    load([slist{i} CommonFileName],'tsIn','tpIn','lapseTrials','B','Fit','WM')
    
    wms(i) = mean(WM(:,Fit.modelUsed));
    
    TS{i} = tsIn{1}(~lapseTrials{1});% - mean(B(:,Fit.modelUsed));
    TP{i} = tpIn{1}(~lapseTrials{1});% - mean(B(:,Fit.modelUsed));
    Theta{i} = regress(TP{i},[ones(size(TS{i})), TS{i}]);
    res{i} = TP{i} - (Theta{i}(2)*TS{i} + Theta{i}(1));
    
    for j = 1:length(tss)
        mRes(i,j) = mean(res{i}(TS{i} == tss(j)));
        stdE(i,j) = std(res{i}(TS{i} == tss(j)))/sqrt(length(res{i}(TS{i} == tss(j))));
        
        % Test if residuals greater than zero
        [~,pvalTemp] = ttest(res{i}(TS{i} == tss(j)),0,'tail','right');
        pvalGreater(i,j) = pvalTemp;
        
        % Test if residuals are less than zero
        [~,pvalTemp] = ttest(res{i}(TS{i} == tss(j)),0,'tail','left');
        pvalLesser(i,j) = pvalTemp;
    end
    
    TS2{i} = tsIn{2}(~lapseTrials{2});% - mean(B(:,Fit.modelUsed));
    TP2{i} = tpIn{2}(~lapseTrials{2});% - mean(B(:,Fit.modelUsed));
    Theta2{i} = regress(TP2{i},[ones(size(TS2{i})), TS2{i}]);
    res2{i} = TP2{i} - (Theta2{i}(2)*TS2{i} + Theta2{i}(1));
    
    for j = 1:length(tss)
        mRes2(i,j) = mean(res2{i}(TS2{i} == tss(j)));
        stdE2(i,j) = std(res2{i}(TS2{i} == tss(j)))/sqrt(length(res2{i}(TS2{i} == tss(j))));
        
        % Test if residuals greater than zero
        [~,pvalTemp] = ttest(res2{i}(TS2{i} == tss(j)),0,'tail','right');
        pvalGreater2(i,j) = pvalTemp;
        
        % Test if residuals are less than zero
        [~,pvalTemp] = ttest(res2{i}(TS2{i} == tss(j)),0,'tail','left');
        pvalLesser2(i,j) = pvalTemp;
    end
end

%% Find expected residual data
tsvec = tss(1):10:tss(end);
for i = 1:length(ExampleSubject)
    ta = ta_expectation3(tsvec',wms(ExampleSubject(i)),1,TAexpectation.dt,'method','numerical',...
        'trials',simulationN,'Support',[min(tss) max(tss)],'wp',0);
    q = regress(ta',[ones(length(tsvec),1),tsvec(:)]);
    expRes(:,i) = ta - (q(2)*tsvec+q(1));
    
    ta2 = ta_expectation3(tsvec',wms(ExampleSubject(i)),2,TAexpectation.dt,'method','numerical',...
        'trials',simulationN,'Support',[min(tss) max(tss)],'wp',0);
    q = regress(ta2',[ones(length(tsvec),1),tsvec(:)]);
    expRes2(:,i) = ta2 - (q(2)*tsvec+q(1));
    
    % taL = ta_expectation3(tsvec',wms(ExampleSubject),1,TAexpectation.dt,'method','numerical',...
    %     'trials',simulationN,'Support',[min(tss) max(tss)],'wp',0,'Type','WeightedLinear');
    % qL = regress(taL',[ones(length(tsvec),1),tsvec(:)]);
    % expResL = taL - (qL(2)*tsvec+qL(1));
end

%% Split into train and test data, fit polynomial
% for i = 1:length(slist)
%     train = false(size(TS{i}));
%     inds = randsample(length(TS{i}),ceil(length(TS{i})/2));
%     train(inds) = true;
%     test = ~train;
%     
%     trainTS = TS{i}(train);
%     trainTP = TP{i}(train);
%     p1 = polyfit(trainTS,trainTP,1);
%     p3 = polyfit(trainTS,trainTP,3);
%     
%     testTS = TS{i}(test);
%     testTP = TS{i}(test);
%     TPhat1 = polyval(p1,testTS);
%     TPhat3 = polyval(p3,testTS);
%     rmse(i,1) = sqrt(mean( (testTP - TPhat1).^2 ));
%     rmse(i,2) = sqrt(mean( (testTP - TPhat3).^2 ));
% end

%% Plotting
switch Plot
    case {'Yes','Y','y','yes',true,1}
        % Set up colors
        colors =[ 0.3010    0.7450    0.9330;...
            0.8500    0.3250    0.0980;...
            0.9290    0.6940    0.1250;...
            0.4940    0.1840    0.5560;...
            0.4660    0.6740    0.1880;...
            0         0.4470    0.7410;...
            0.6350    0.0780    0.1840;...
            0         0         1.0000;...
            0         0.5000    0;...
            1.0000    0         0];
        
        figure('Name','Mean residuals','Position',[82 533 1748 333])
        subplot(1,2,1)
        for i = 1:length(slist)
            %plot(TS{i},res{i},'o','Color',[0.7 0.7 0.7]);
            hold on
            if ~any(i == ExampleSubject)
                errorbar(tss,mRes(i,:),stdE(i,:),'.','LineWidth',2,'Color',colors(i,:));
                h(i) = plot(tss,mRes(i,:),'o','Color','none');
                set(h(i),'MarkerFaceColor',colors(i,:));
            end
        end
        for i = 1:length(ExampleSubject)
            errorbar(tss,mRes(ExampleSubject(i),:),stdE(ExampleSubject(i),:),...
                '.','LineWidth',2,'Color',colors(ExampleSubject(i),:));
            h(ExampleSubject(i)) = plot(tss,mRes(ExampleSubject(i),:),'o','Color','none');
            set(h(ExampleSubject(i)),'MarkerFaceColor',colors(ExampleSubject(i),:));
            plot(tsvec,expRes(:,i),'-','LineWidth',2,'Color',colors(ExampleSubject(i),:))
        end
        %plot(tsvec,expResL,'--','LineWidth',2,'Color',colors(ExampleSubject,:))
        plotHorizontal(0);
%         ax = axis;
%         for i = 1:length(slist)
%             for j = 1:length(tss)
%                 if pvalGreater(i,j) <= alpha
%                     plot(tss(j),ax(4)+i-1,'*','Color',colors(i,:));
%                 end
%                 if pvalLesser(i,j) <= alpha
%                     plot(tss(j),ax(3)-i-1,'*','Color',colors(i,:));
%                 end
%             end
%         end     
        xlabel('t_s')
        ylabel('Residual')
        legend(h,slist)
        mymakeaxis(gca,'xytitle','RSG','xticks',[600 800 1000],...
            'xticklabels',{'600','800','1000'})
        
        
        subplot(1,2,2)
        for i = 1:length(slist)
            %plot(TS{i},res{i},'o','Color',[0.7 0.7 0.7]);
            hold on
            if ~any(i == ExampleSubject)
                errorbar(tss,mRes2(i,:),stdE2(i,:),'.','LineWidth',2,'Color',colors(i,:));
                h(i) = plot(tss,mRes2(i,:),'o','Color','none');
                set(h(i),'MarkerFaceColor',colors(i,:));
            end
        end
        for i = 1:length(ExampleSubject)
            errorbar(tss,mRes2(ExampleSubject(i),:),stdE2(ExampleSubject(i),:),...
                '.','LineWidth',2,'Color',colors(ExampleSubject(i),:));
            h(ExampleSubject(i)) = plot(tss,mRes2(ExampleSubject(i),:),'o','Color','none');
            set(h(ExampleSubject(i)),'MarkerFaceColor',colors(ExampleSubject(i),:));
            plot(tsvec,expRes2(:,i),'-','LineWidth',2,'Color',colors(ExampleSubject(i),:))
        end
        %plot(tsvec,expResL,'--','LineWidth',2,'Color',colors(ExampleSubject,:))
        plotHorizontal(0);
%         ax = axis;
%         for i = 1:length(slist)
%             for j = 1:length(tss)
%                 if pvalGreater(i,j) <= alpha
%                     plot(tss(j),ax(4)+i-1,'*','Color',colors(i,:));
%                 end
%                 if pvalLesser(i,j) <= alpha
%                     plot(tss(j),ax(3)-i-1,'*','Color',colors(i,:));
%                 end
%             end
%         end     
        xlabel('t_s')
        ylabel('Residual')
        legend(h,slist)
        mymakeaxis(gca,'xytitle','RSSG','xticks',[600 800 1000],...
            'xticklabels',{'600','800','1000'})
        
    case {'No','no','N','n',false,0}
        
    otherwise
        error(['Plot option ' Plot ' not recognized!'])
        
end