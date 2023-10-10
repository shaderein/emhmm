function [px, Xgrid, Ygrid] = gmm2heatmap(gmm, Xp, Yp)
% gmm2heatmap - convert a GMM to a heatmap representation
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-21
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% VERSIONS
% 2021-07-21: v0.78 - ABC - initial version


% collect all points
[Xgrid,Ygrid] = meshgrid(Xp, Yp);

Xv = Xgrid(:);
Yv = Ygrid(:);

y = [Xv, Yv];

% compute probabilities
%P = zeros(size(y,1),1);
%
%for k=1:length(gmm.pdf)
%    myp = mvnpdf(y, gmm.pdf{k}.mean, gmm.pdf{k}.cov);
%    P = P + gmm.prior(k)*myp;
%end

P = gmm_pdf(gmm, y);

% make sure it sums to 1
P = P/sum(P(:));

% reshape
px = reshape(P, size(Xgrid));

return

%% test
pdfy = vbhmm_prob_marginalobs(hmms{4}, 6)
figure(100)
subplot(2,4,1)
imagesc(img);
axis image
subplot(2,4,2)
plot_emissions([], [], pdfy{1}.pdf, img)
for t=1:length(pdfy)
    subplot(2,4,t+2)
    [py, Xg, Yg] = gmm2heatmap(pdfy{t}, [1 512], [1 384], [512 384]);
    imagesc(Xg(1, [1 end]), Yg([1 end], 1), py)
    axis image
end
