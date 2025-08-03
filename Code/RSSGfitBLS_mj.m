function d = RSSGfitBLS_mj(d)
 
%%
% fit according to BLS
d.OPTIONS = optimset('Display','iter');
BLSminimizant = @(p) get_pRSSG(p(1), p(2), d);
 
d.minimizer = 'fminsearch(BLSminimizant, [wm2fit_ini wintfit_ini], d.OPTIONS);';
wm2fit_ini = d.wm1;
wintfit_ini = d.wm1;
[d.lparams, d.llikelihood] = eval(d.minimizer);
 
%%
% save results
d.wm2RSSG_fit = d.lparams(1);
d.wintRSSG_fit = d.lparams(2);
%save(d.filepath, '-struct', 'd');
 
 
function p_tp_ts = get_pRSSG(wm2fit, wintfit, d)
%%
tm1 = d.tmin-5*d.wm1*d.tmin:ceil(d.wm1*20):d.tmax+5*d.wm1*d.tmax;
tm2 = d.tmin-5*wm2fit*d.tmin:ceil(wm2fit*20):d.tmax+5*wm2fit*d.tmax;
[TM1, TM2] = meshgrid(tm1,tm2);
 
[~, TM_MLE] = RSSG_MLE(TM1, TM2, d.wm1, wm2fit);
fBLS = sim_RSG_RSSG_fBLS(d.tmin, d.tmax, wintfit, 0, 'RSG');
TE = fBLS(TM_MLE(:)');
TE = reshape(TE',size(TM_MLE));
Lhood = zeros(1,length(d.ts));
for i = 1:length(d.ts)
    prob1 = exp(-1./2.*(d.ts(i)-TM1).^2./d.wm1.^2./d.ts(i).^2);
    prob2 = exp(-1./2.*(d.ts(i)-TM2).^2./wm2fit.^2./d.ts(i).^2);
    prob_p = exp(-1./2.*(d.tp(i)-TE).^2./d.wp.^2./TE.^2)./TE;
    F = prob1.*prob2.*prob_p;
    Lhood(i) = trapz(tm2,trapz(tm1,F,2))./d.ts(i)./d.ts(i)./wm2fit;
end
p_tp_ts = -sum(log(Lhood));
 
function fBLS = sim_RSG_RSSG_fBLS(tmin, tmax, wM, K, condition)
 
%%
tsvec = tmin:tmax;
%%
switch condition
    case 'RSG'
        tm1vec = tmin-7*wM*tmin:ceil(wM*20):tmax+7*wM*tmax;
        [xx, zz] = meshgrid(tm1vec,tsvec);
        v1 = xx-zz;
        s11 = 1./wM.^2./zz.^2;
        I1 = exp(-1./2.*s11.*v1.^2);
        I2 = exp(-1./2.*s11.*v1.^2)./zz;
        f_tmvec = trapz(I1,1)./trapz(I2,1);
        fBLS = @(tm) interp1(tm1vec,f_tmvec,tm);
    case 'RSSG'
        tm1vec = tmin-3*wM*tmin:tmax+3*wM*tmax;
        tm2vec = tmin-3*wM*tmin:tmax+3*wM*tmax;
        [xx, yy, zz] = meshgrid(tm1vec,tm2vec,tsvec);
        v1 = xx-zz;
        v2 = yy-zz;
        s11 = 1./wM.^2./zz.^2;
        I1 = exp(-1./2.*s11.*(v1.^2+v2.^2+2*K*v1.*v2))./zz;
        I2 = exp(-1./2.*s11.*(v1.^2+v2.^2+2*K*v1.*v2))./zz.^2;
        f_tmvec = trapz(I1,3)./trapz(I2,3);
        fBLS = @(tm1,tm2) interp2(tm1vec,tm2vec,f_tmvec,tm1,tm2);
end
 
function [w12, t_mle] = RSSG_MLE(tm1, tm2, w1, w2)
 
a1 = w2^2/(w1.^2+w2.^2);
a2 = w1^2/(w1.^2+w2.^2);
w12 = w1*w2/sqrt(w1.^2+w2.^2);
mhat = (a1*tm1+a2*tm2)/2;
m2hat = (a1*tm1.^2+a2*tm2.^2)/2;
t_mle = mhat.*(-1+sqrt(1+4*w12.^2.*m2hat./mhat.^2))./2./w12.^2;