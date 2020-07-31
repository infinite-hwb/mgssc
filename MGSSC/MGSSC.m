clc
clear all
addpath(genpath(cd))
addpath('.\dataset');
addpath('.\performance');
addpath('.\tool');
load ORL
truth = gnd;
k = max(gnd);
K = KERNEL;
n = size(fea,1);
%% MGSSC
fprintf('\nCentroid multiview MGSSC\n');
for mu=[2]
    for lambda1=[0.5]
        for lambda2=[0.5]
            lambda = [lambda1 lambda2];
            [A,C] = centroid_MGSSC(n,K,mu,lambda); % joint affinity matrix
            fprintf('mu is %g \n' , mu) ;
            fprintf('lambda1 is %g \n' , lambda1) ;
            fprintf('labmda2 is %g \n' , lambda2) ;
            [best.CA best.F best.P best.R best.nmi best.AR] = spectral_clustering(A, k, truth);
            best
        end
    end
end