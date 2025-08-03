function Theta = printParametersFullLLvLNOCV(slist,varargin)
%% TestLinearity
%
%   h = TestLinearity(slist)
%
%   Tests the hypothesis that reproduction data is linearly related to
%   sample interval.
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'slist')     % List of subjects to analyze
addParameter(Parser,'Models',{'LNE','EKF','BLS','BLS_wm1wm2'})

parse(Parser,slist,varargin{:})

slist = Parser.Results.slist;
Models = Parser.Results.Models;

Theta.models = Models;
Theta.subjects = slist;
%% Sweep through modles/subjects and record the information
for modeli = 1:length(Models)
    for si = 1:length(slist)
        switch Models{modeli}
            
            case 'LNE'
                % LNOCV
                clear WM WP B lapse WM_DRIFT
                load([slist{si} '_BLSbiasedFitResults20150913'],'WM','WP','B','lapse','bias','BIASs',...
                    'variance','VARs','rmse','RMSEs','Fit','Llikelihood','mtp_in',...
                    'notFitLlikelihood','stdtp_in','WM_DRIFT','W_INT','ALPHA',...
                    'tsIn','lapseTrials')
                
                modelInd = strcmp(Fit.fittype,'aveMeasurements');
                % Parameters fit to data
                Theta.LNOCV.N(si,:,modeli) = size(WM(:,modelInd),1);
                Theta.LNOCV.wm(si,:,modeli) = mean(WM(:,modelInd),1);
                Theta.LNOCV.wp(si,:,modeli) = mean(WP(:,modelInd,1),1);
                Theta.LNOCV.b(si,:,modeli) = mean(B(:,modelInd),1);
                Theta.LNOCV.Lapse(si,:,modeli) = mean(lapse(:,modelInd),1);
                Theta.LNOCV.LL = -Llikelihood(1:end-1,modelInd);
                if exist('WM_DRIFT','var')
                    if size(WM_DRIFT,2) == 1
                        Theta.LNOCV.wm_drift(si,:,modeli) = [mean(WM_DRIFT(:,modelInd),1) NaN];
                    else
                        Theta.LNOCV.wm_drift(si,:,modeli) = mean(WM_DRIFT(:,modelInd),1);
                    end
                else
                    Theta.LNOCV.wm_drift(si,:,modeli) = [NaN NaN];
                end
                if exist('W_INT','var')
                    if size(W_INT,2) == 1
                        Theta.LNOCV.w_int(si,:,modeli) = [mean(W_INT(:,modelInd),1) NaN];
                    else
                        Theta.LNOCV.w_int(si,:,modeli) = mean(W_INT(:,modelInd),1);
                    end
                else
                    Theta.LNOCV.w_int(si,:,modeli) = [NaN NaN];
                end
                if exist('ALPHA','var')
                    if size(ALPHA,2) == 1
                        Theta.LNOCV.alpha(si,:,modeli) = [mean(ALPHA(:,modelInd),1) NaN];
                    else
                        Theta.LNOCV.alpha(si,:,modeli) = mean(ALPHA(:,modelInd),1);
                    end
                else
                    Theta.LNOCV.alpha(si,:,modeli) = [NaN NaN];
                end
                
                
                clear WM WP B lapse WM_DRIFT                
                % Full LL
                load([slist{si} '_BLSbiasedLapse_ObsAct0_20180625'],'WM','WP','B','lapse','bias','BIASs',...
                    'variance','VARs','rmse','RMSEs','Fit','Llikelihood','mtp_in',...
                    'notFitLlikelihood','stdtp_in','WM_DRIFT','W_INT','ALPHA',...
                    'tsIn','lapseTrials')
                
                modelInd = strcmp(Fit.fittype,'aveMeasurements');
                % Parameters fit to data
                Theta.Full.N(si,:,modeli) = size(WM(:,modelInd),1);
                Theta.Full.wm(si,:,modeli) = WM(:,modelInd);
                Theta.Full.wp(si,:,modeli) = WP(:,modelInd,1);
                Theta.Full.b(si,:,modeli) = B(:,modelInd);
                Theta.Full.Lapse(si,:,modeli) = lapse(:,modelInd);
                Theta.Full.LL = -Llikelihood(1:end-1,modelInd);
                
                if exist('WM_DRIFT','var')
                    if size(WM_DRIFT,2) == 1
                        Theta.Full.wm_drift(si,:,modeli) = [WM_DRIFT(:,modelInd) NaN];
                    else
                        Theta.Full.wm_drift(si,:,modeli) = WM_DRIFT(:,modelInd);
                    end
                else
                    Theta.Full.wm_drift(si,:,modeli) = [NaN NaN];
                end
                if exist('W_INT','var')
                    if size(W_INT,2) == 1
                        Theta.Full.w_int(si,:,modeli) = [W_INT(:,modelInd) NaN];
                    else
                        Theta.Full.w_int(si,:,modeli) = W_INT(:,modelInd);
                    end
                else
                    Theta.Full.w_int(si,:,modeli) = [NaN NaN];
                end
                if exist('ALPHA','var')
                    if size(ALPHA,2) == 1
                        Theta.Full.alpha(si,:,modeli) = [ALPHA(:,modelInd) NaN];
                    else
                        Theta.Full.alpha(si,:,modeli) = ALPHA(:,modelInd);
                    end
                else
                    Theta.Full.alpha(si,:,modeli) = [NaN NaN];
                end
                
            case 'EKF'
                % LNOCV
                clear WM WP B lapse WM_DRIFT
                load([slist{si} '_EKF_ObsAct0_20171125'],'WM','WP','B','lapse','bias','BIASs',...
                    'variance','VARs','rmse','RMSEs','Fit','Llikelihood','mtp_in',...
                    'notFitLlikelihood','stdtp_in','WM_DRIFT','W_INT','ALPHA',...
                    'tsIn','lapseTrials')
                
                modelInd = strcmp(Fit.fittype,'EKF');
                % Parameters fit to data
                Theta.LNOCV.N(si,:,modeli) = size(WM(:,modelInd),1);
                Theta.LNOCV.wm(si,:,modeli) = mean(WM(:,modelInd),1);
                Theta.LNOCV.wp(si,:,modeli) = mean(WP(:,modelInd,1),1);
                Theta.LNOCV.b(si,:,modeli) = mean(B(:,modelInd),1);
                Theta.LNOCV.Lapse(si,:,modeli) = mean(lapse(:,modelInd),1);
                Theta.LNOCV.LL = -Llikelihood(1:end-1,modelInd);
                if exist('WM_DRIFT','var')
                    if size(WM_DRIFT,2) == 1
                        Theta.LNOCV.wm_drift(si,:,modeli) = [mean(WM_DRIFT(:,modelInd),1) NaN];
                    else
                        Theta.LNOCV.wm_drift(si,:,modeli) = mean(WM_DRIFT(:,modelInd),1);
                    end
                else
                    Theta.LNOCV.wm_drift(si,:,modeli) = [NaN NaN];
                end
                if exist('W_INT','var')
                    if size(W_INT,2) == 1
                        Theta.LNOCV.w_int(si,:,modeli) = [mean(W_INT(:,modelInd),1) NaN];
                    else
                        Theta.LNOCV.w_int(si,:,modeli) = mean(W_INT(:,modelInd),1);
                    end
                else
                    Theta.LNOCV.w_int(si,:,modeli) = [NaN NaN];
                end
                if exist('ALPHA','var')
                    if size(ALPHA,2) == 1
                        Theta.LNOCV.alpha(si,:,modeli) = [mean(ALPHA(:,modelInd),1) NaN];
                    else
                        Theta.LNOCV.alpha(si,:,modeli) = mean(ALPHA(:,modelInd),1);
                    end
                else
                    Theta.LNOCV.alpha(si,:,modeli) = [NaN NaN];
                end
                
                
                clear WM WP B lapse WM_DRIFT                
                % Full LL
                load([slist{si} '_EKF_ObsAct0_20180625'],'WM','WP','B','lapse','bias','BIASs',...
                    'variance','VARs','rmse','RMSEs','Fit','Llikelihood','mtp_in',...
                    'notFitLlikelihood','stdtp_in','WM_DRIFT','W_INT','ALPHA',...
                    'tsIn','lapseTrials')
                
                modelInd = strcmp(Fit.fittype,'EKF');
                % Parameters fit to data
                Theta.Full.N(si,:,modeli) = size(WM(:,modelInd),1);
                Theta.Full.wm(si,:,modeli) = WM(:,modelInd);
                Theta.Full.wp(si,:,modeli) = WP(:,modelInd,1);
                Theta.Full.b(si,:,modeli) = B(:,modelInd);
                Theta.Full.Lapse(si,:,modeli) = lapse(:,modelInd);
                Theta.Full.LL = -Llikelihood(1:end-1,modelInd);
                
                if exist('WM_DRIFT','var')
                    if size(WM_DRIFT,2) == 1
                        Theta.Full.wm_drift(si,:,modeli) = [WM_DRIFT(:,modelInd) NaN];
                    else
                        Theta.Full.wm_drift(si,:,modeli) = WM_DRIFT(:,modelInd);
                    end
                else
                    Theta.Full.wm_drift(si,:,modeli) = [NaN NaN];
                end
                if exist('W_INT','var')
                    if size(W_INT,2) == 1
                        Theta.Full.w_int(si,:,modeli) = [W_INT(:,modelInd) NaN];
                    else
                        Theta.Full.w_int(si,:,modeli) = W_INT(:,modelInd);
                    end
                else
                    Theta.Full.w_int(si,:,modeli) = [NaN NaN];
                end
                if exist('ALPHA','var')
                    if size(ALPHA,2) == 1
                        Theta.Full.alpha(si,:,modeli) = [ALPHA(:,modelInd) NaN];
                    else
                        Theta.Full.alpha(si,:,modeli) = ALPHA(:,modelInd);
                    end
                else
                    Theta.Full.alpha(si,:,modeli) = [NaN NaN];
                end
                
                
            case 'BLS'
                % LNOCV
                clear WM WP B lapse WM_DRIFT
                load([slist{si} '_BLSbiasedFitResults20150913'],'WM','WP','B','lapse','bias','BIASs',...
                    'variance','VARs','rmse','RMSEs','Fit','Llikelihood','mtp_in',...
                    'notFitLlikelihood','stdtp_in','WM_DRIFT','W_INT','ALPHA',...
                    'tsIn','lapseTrials')
                
                modelInd = strcmp(Fit.fittype,'BLSbiasedLapse');
                % Parameters fit to data
                Theta.LNOCV.N(si,:,modeli) = size(WM(:,modelInd),1);
                Theta.LNOCV.wm(si,:,modeli) = mean(WM(:,modelInd),1);
                Theta.LNOCV.wp(si,:,modeli) = mean(WP(:,modelInd,1),1);
                Theta.LNOCV.b(si,:,modeli) = mean(B(:,modelInd),1);
                Theta.LNOCV.Lapse(si,:,modeli) = mean(lapse(:,modelInd),1);
                Theta.LNOCV.LL = -Llikelihood(1:end-1,modelInd);
                if exist('WM_DRIFT','var')
                    if size(WM_DRIFT,2) == 1
                        Theta.LNOCV.wm_drift(si,:,modeli) = [mean(WM_DRIFT(:,modelInd),1) NaN];
                    else
                        Theta.LNOCV.wm_drift(si,:,modeli) = mean(WM_DRIFT(:,modelInd),1);
                    end
                else
                    Theta.LNOCV.wm_drift(si,:,modeli) = [NaN NaN];
                end
                if exist('W_INT','var')
                    if size(W_INT,2) == 1
                        Theta.LNOCV.w_int(si,:,modeli) = [mean(W_INT(:,modelInd),1) NaN];
                    else
                        Theta.LNOCV.w_int(si,:,modeli) = mean(W_INT(:,modelInd),1);
                    end
                else
                    Theta.LNOCV.w_int(si,:,modeli) = [NaN NaN];
                end
                if exist('ALPHA','var')
                    if size(ALPHA,2) == 1
                        Theta.LNOCV.alpha(si,:,modeli) = [mean(ALPHA(:,modelInd),1) NaN];
                    else
                        Theta.LNOCV.alpha(si,:,modeli) = mean(ALPHA(:,modelInd),1);
                    end
                else
                    Theta.LNOCV.alpha(si,:,modeli) = [NaN NaN];
                end
                
                
                clear WM WP B lapse WM_DRIFT                
                % Full LL
                load([slist{si} '_BLSbiasedLapse_ObsAct0_20180625'],'WM','WP','B','lapse','bias','BIASs',...
                    'variance','VARs','rmse','RMSEs','Fit','Llikelihood','mtp_in',...
                    'notFitLlikelihood','stdtp_in','WM_DRIFT','W_INT','ALPHA',...
                    'tsIn','lapseTrials')
                
                modelInd = strcmp(Fit.fittype,'BLSbiasedLapse');
                % Parameters fit to data
                Theta.Full.N(si,:,modeli) = size(WM(:,modelInd),1);
                Theta.Full.wm(si,:,modeli) = WM(:,modelInd);
                Theta.Full.wp(si,:,modeli) = WP(:,modelInd,1);
                Theta.Full.b(si,:,modeli) = B(:,modelInd);
                Theta.Full.Lapse(si,:,modeli) = lapse(:,modelInd);
                Theta.Full.LL = -Llikelihood(1:end-1,modelInd);
                
                if exist('WM_DRIFT','var')
                    if size(WM_DRIFT,2) == 1
                        Theta.Full.wm_drift(si,:,modeli) = [WM_DRIFT(:,modelInd) NaN];
                    else
                        Theta.Full.wm_drift(si,:,modeli) = WM_DRIFT(:,modelInd);
                    end
                else
                    Theta.Full.wm_drift(si,:,modeli) = [NaN NaN];
                end
                if exist('W_INT','var')
                    if size(W_INT,2) == 1
                        Theta.Full.w_int(si,:,modeli) = [W_INT(:,modelInd) NaN];
                    else
                        Theta.Full.w_int(si,:,modeli) = W_INT(:,modelInd);
                    end
                else
                    Theta.Full.w_int(si,:,modeli) = [NaN NaN];
                end
                if exist('ALPHA','var')
                    if size(ALPHA,2) == 1
                        Theta.Full.alpha(si,:,modeli) = [ALPHA(:,modelInd) NaN];
                    else
                        Theta.Full.alpha(si,:,modeli) = ALPHA(:,modelInd);
                    end
                else
                    Theta.Full.alpha(si,:,modeli) = [NaN NaN];
                end
                
                
                
            case 'BLS_wm1wm2'
                % LNOCV
                clear WM WP B lapse WM_DRIFT
                load([slist{si} '_BLS_wm1wm2_ObsAct0_20171024'],'WM','WP','B','lapse','bias','BIASs',...
                    'variance','VARs','rmse','RMSEs','Fit','Llikelihood','mtp_in',...
                    'notFitLlikelihood','stdtp_in','WM_DRIFT','W_INT','ALPHA',...
                    'tsIn','lapseTrials')
                
                modelInd = strcmp(Fit.fittype,'BLS_wm1wm2');
                % Parameters fit to data
                Theta.LNOCV.N(si,:,modeli) = size(WM(:,modelInd),1);
                Theta.LNOCV.wm(si,:,modeli) = mean(WM(:,modelInd),1);
                Theta.LNOCV.wp(si,:,modeli) = mean(WP(:,modelInd,1),1);
                Theta.LNOCV.b(si,:,modeli) = mean(B(:,modelInd),1);
                Theta.LNOCV.Lapse(si,:,modeli) = mean(lapse(:,modelInd),1);
                Theta.LNOCV.LL = -Llikelihood(1:end-1,modelInd);
                if exist('WM_DRIFT','var')
                    if size(WM_DRIFT,2) == 1
                        Theta.LNOCV.wm_drift(si,:,modeli) = [mean(WM_DRIFT(:,modelInd),1) NaN];
                    else
                        Theta.LNOCV.wm_drift(si,:,modeli) = mean(WM_DRIFT(:,modelInd),1);
                    end
                else
                    Theta.LNOCV.wm_drift(si,:,modeli) = [NaN NaN];
                end
                if exist('W_INT','var')
                    if size(W_INT,2) == 1
                        Theta.LNOCV.w_int(si,:,modeli) = [mean(W_INT(:,modelInd),1) NaN];
                    else
                        Theta.LNOCV.w_int(si,:,modeli) = mean(W_INT(:,modelInd),1);
                    end
                else
                    Theta.LNOCV.w_int(si,:,modeli) = [NaN NaN];
                end
                if exist('ALPHA','var')
                    if size(ALPHA,2) == 1
                        Theta.LNOCV.alpha(si,:,modeli) = [mean(ALPHA(:,modelInd),1) NaN];
                    else
                        Theta.LNOCV.alpha(si,:,modeli) = mean(ALPHA(:,modelInd),1);
                    end
                else
                    Theta.LNOCV.alpha(si,:,modeli) = [NaN NaN];
                end
                
                                
                % Full LL
                clear WM WP B lapse WM_DRIFT
                load([slist{si} '_BLS_wm1wm2_ObsAct0_20180625'],'WM','WP','B','lapse','bias','BIASs',...
                    'variance','VARs','rmse','RMSEs','Fit','Llikelihood','mtp_in',...
                    'notFitLlikelihood','stdtp_in','WM_DRIFT','W_INT','ALPHA',...
                    'tsIn','lapseTrials')
                
                modelInd = strcmp(Fit.fittype,'BLS_wm1wm2');
                % Parameters fit to data
                Theta.Full.N(si,:,modeli) = size(WM(:,modelInd),1);
                Theta.Full.wm(si,:,modeli) = WM(:,modelInd);
                Theta.Full.wp(si,:,modeli) = WP(:,modelInd,1);
                Theta.Full.b(si,:,modeli) = B(:,modelInd);
                Theta.Full.Lapse(si,:,modeli) = lapse(:,modelInd);
                Theta.Full.LL = -Llikelihood(1:end-1,modelInd);
                
                if exist('WM_DRIFT','var')
                    if size(WM_DRIFT,2) == 1
                        Theta.Full.wm_drift(si,:,modeli) = [WM_DRIFT(:,modelInd) NaN];
                    else
                        Theta.Full.wm_drift(si,:,modeli) = WM_DRIFT(:,modelInd);
                    end
                else
                    Theta.Full.wm_drift(si,:,modeli) = [NaN NaN];
                end
                if exist('W_INT','var')
                    if size(W_INT,2) == 1
                        Theta.Full.w_int(si,:,modeli) = [W_INT(:,modelInd) NaN];
                    else
                        Theta.Full.w_int(si,:,modeli) = W_INT(:,modelInd);
                    end
                else
                    Theta.Full.w_int(si,:,modeli) = [NaN NaN];
                end
                if exist('ALPHA','var')
                    if size(ALPHA,2) == 1
                        Theta.Full.alpha(si,:,modeli) = [ALPHA(:,modelInd) NaN];
                    else
                        Theta.Full.alpha(si,:,modeli) = ALPHA(:,modelInd);
                    end
                else
                    Theta.Full.alpha(si,:,modeli) = [NaN NaN];
                end
        end
        
    end  % end subject loop
end % end model loop