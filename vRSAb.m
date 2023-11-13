function vbRSA_boot(directroot,betaroot,designroot, contrastroot,outroot)
%last update - 11-28-22 7:40am - STG (resized b_dist)


%% user inputs 
num_boots = 1000;
%directroot = directory identifier, in example I use 'bd2';
%betaroot= beta file identifier, in example I use 'Precentral_L';
%designroot= design matrix identifier, in example I use 'design';
%contrastroot= identifier for file containing what contrasts to do, in example I use 'contrast'
%outroot=something in the output files to identify them

%Ns=4;

%%Scott's loader code, amended to fit bootstrap conventions%%
clear data design
subs=dir(strcat(directroot,'_*'));
[Nn s]=size(subs);
nan_ind=0;
for sub=1:Nn
    study = dir(fullfile(subs(sub).folder, subs(sub).name, strcat(directroot,'*_',betaroot,'.txt')));
    Nr(sub)=size(study,1);
    sub_run_data=[];
    for run=1:Nr(sub)
        sub_run = load(fullfile(study(run).folder, study(run).name ),'-ascii' );
        betas_per_run(run,sub)=size(sub_run,1)-3; %index needed to deal with missing betas
        sub_run_data = cat(1,sub_run_data,sub_run(4:end,:));
    end
    data{sub}=sub_run_data;
    nan_ind=nan_ind+sum(isnan(sub_run_data),1);
    
    %and load design matrix
    des1 = dir(fullfile(subs(sub).folder, subs(sub).name, strcat(designroot,'*.txt')));
    Xtemp = load( fullfile(des1.folder, des1.name ),'-ascii' );
    design{sub}=Xtemp;

    % and make nuisance regressors (related to run)
    nuisance{sub}=dummyvar(repelem(1:Nr(sub),betas_per_run(:,sub)));
end

%remove voxels if any subject had a NaN there
keep_vox=nan_ind==0;
voxels_kept=sum(keep_vox);
disp(['keeping ',num2str(voxels_kept),' voxels']);
data=cellfun(@(x) x(:,keep_vox),data,'UniformOutput',false);

Nvox = size(data{1},2);

% will look for a mat file called "contrasts"
% f-by-c array
% a row for each factor (column in the design matrix) (f)
% a column for each contrast (c)
contrasts=load(strcat(contrastroot,'.txt'),'-ascii' );
%contrasts=load('contrast.txt');

%Nn = size(data,2);
%Np = size(data{1},2);

ortho_contrasts      = spm_orth(contrasts,'norm');       % orthonormalise
Nc     = size(ortho_contrasts,2);
for i = 1:Nc
    C{i} = ortho_contrasts(:,i)*ortho_contrasts(:,i)';
end

%Prep for the bootstrap
all_data_2d = [];
%max_rows = cellfun(@(x) size(x,1),data);
%padding = sum(max(max_rows)-max_rows);
for sub = 1:Nn;
    all_data_2d = cat(1,all_data_2d,data{sub});
end

NnNm = size(all_data_2d,1);
b_dist = zeros(num_boots,size(ortho_contrasts,2)+1);
M.X       = ones(Nn,1);
%% compute null dist
tic;
for b = 1:num_boots
    
    clear RSA
    for sub = 1:Nn;

        %Extract design and nuisance matrix from cell arrays
        Z = design{sub};
        B = nuisance{sub};
        
        %Permute a data matrix, same number of measuremebts
        y = all_data_2d(randperm(NnNm,size(Z,1)),:);
        y = y-mean(y,2);

        %Add additional contrast to new contrast matrix (Q)
        Q = C;
        Q{end+1}=pinv(Z)*pinv(Z)';

        %Compute spatial covariance param (Nv)
        R = eye(size(Z,1)) - [Z B]*pinv([Z B]);
        e  = R*y;
        e  = e'*e;
        Nv = trace(e)^2/trace(e*e);

        %compute matrix for confounds
        conf_mat = pinv(Z)*B;

        %finally pattern covariance matrix G = U*U'
        U = pinv(Z)*y;
        G = U*U';

        [Cy,h,Ph,F,Fa,Fc,Eh,Ch,hE,hC,Qh] = spm_reml_sc(G,conf_mat,Q,Nv,-16,128);
        % accumulate these for input into next step
        %----------------------------------------------------------------------
        RSA{sub}.M.pE = hE;         % prior expectation of parameters
        RSA{sub}.M.pC = hC;         % prior covariances of parameters
        RSA{sub}.Ep   = Eh;         % posterior expectations
        RSA{sub}.Cp   = Ch;         % posterior covariance
        RSA{sub}.Q    = Qh;         % scaled covariance components
        RSA{sub}.F    = F;          % free energy

    end

    % parametric empirical Bayes to compute group average hyper parameters
    [PEB,RSA] = spm_dcm_peb(RSA(:),M);

    % Bayesian model reduction for model comparison
    %==========================================================================
    pE    = PEB.M.pE;
    pC    = PEB.M.pC;
    qE    = PEB.Ep;
    qC    = PEB.Cp;

    for i = 1:Nc+1
        % Place precise shrinkage priors on each component
        %----------------------------------------------------------------------
        rC      = pC;
        %we could use a value of 1/128 etc...
        rC(i,i) = 1/128;
        Fc(i,1) = spm_log_evidence(qE,qC,pE,pC,pE,rC);

        % and assess evidence for just this component
        %----------------------------------------------------------------------
        %         rC      = pC;
        %         jj       = 1:Nc+1; jj(i) = [];
        %         rC(jj,jj) = 1/128;
        %         Fs(i,1)    = spm_log_evidence(qE,qC,pE,pC,pE,rC);
    end

    b_dist(b,:) = Fc';

end
toc;



%% now do real data
clear RSA
for sub = 1:Nn;

    %Extract data and design matrix from cell arrays
    y = data{sub};
    y = y-mean(y,2);
    Z = design{sub};
    B = nuisance{sub};

    %Add additional contrast to new contrast matrix (Q)
    Q = C;
    Q{end+1}=pinv(Z)*pinv(Z)';

    %Compute spatial covariance param (Nv)
    R = eye(size(Z,1)) - [Z B]*pinv([Z B]);
    e  = R*y;
    e  = e'*e;
    Nv = trace(e)^2/trace(e*e);

    %compute matrix for confounds (it'll just be a cond-by-2 ones)
    conf_mat = pinv(Z)*B;

    %finally pattern covariance matrix G = U*U'
    U = pinv(Z)*y;
    G = U*U';

    [Cy,h,Ph,F,Fa,Fc,Eh,Ch,hE,hC,Qh] = spm_reml_sc(G,conf_mat,Q,Nv,-16,128);
    % accumulate these for input into next step
    %----------------------------------------------------------------------
    RSA{sub}.M.pE = hE;         % prior expectation of parameters
    RSA{sub}.M.pC = hC;         % prior covariances of parameters
    RSA{sub}.Ep   = Eh;         % posterior expectations
    RSA{sub}.Cp   = Ch;         % posterior covariance
    RSA{sub}.Q    = Qh;         % scaled covariance components
    RSA{sub}.F    = F;          % free energy
end

% parametric empirical Bayes to compute group average hyper parameters
M.X       = ones(Nn,1);
[PEB,RSA] = spm_dcm_peb(RSA(:),M);

% Bayesian model reduction for model comparison
%==========================================================================
pE    = PEB.M.pE;
pC    = PEB.M.pC;
qE    = PEB.Ep;
qC    = PEB.Cp;

for i = 1:Nc+1
    % Place precise shrinkage priors on each component
    %----------------------------------------------------------------------
    rC      = pC;
    %we could use a value of 1/128 etc...
    rC(i,i) = 1/128;
    Fc(i,1) = spm_log_evidence(qE,qC,pE,pC,pE,rC);

    % and assess evidence for just this component
    %----------------------------------------------------------------------
    rC      = pC;
    jj       = 1:Nc+1; jj(i) = [];
    rC(jj,jj) = 1/128;
    Fs(i,1)    = spm_log_evidence(qE,qC,pE,pC,pE,rC);
end

figure('Units','normalized','OuterPosition',[0 0 1 1]);

clf;
subplot(3,6,1);
imagesc(data{1});
ylabel('measurements');
xlabel('voxels');
title('example subject data');
subplot(3,6,7);
imagesc(design{1});
title('design matrix');
ylabel('measurements');
xlabel('conditions');
subplot(3,6,13);
imagesc(nuisance{1});
title('nuisance matrix');
ylabel('measurements');
xlabel('run');

for c = 1:size(Q,2);
    subplot(size(Q,2),4,-2 + 4*c);
    imagesc(Q{c});axis square;title(['contrast #',num2str(c)]);

    subplot(size(Q,2),4,-1 + 4*c);
    histogram(b_dist(:,c),50,'EdgeColor','None');hold on;
    yl = ylim;
    line([Fc(c) Fc(c)],yl,'LineWidth',3,'color','k');
    line([median(b_dist(:,c)) median(b_dist(:,c))],yl,'LineWidth',3,'color',[.5,.5,.5],'LineStyle','--');
    ylim([0 yl(2)*1.2]);
    text(median(b_dist(:,c)),yl(2)*1.1,'Null Dist','HorizontalAlignment','center');
    xlabel('log evidence');
    set(gca,'YTick',[]);
    axis square

    subplot(size(Q,2),4,4*c);
    BF = b_dist(:,c) - Fc(c);
    pts=linspace(0.1,99.9-95,20);
    pt1=prctile(BF,pts);
    pt2=prctile(BF,95+pts);
    cis=abs(pt2-pt1);
    [foo,hpdpi]=min(cis);
    H=[pt1(hpdpi) pt2(hpdpi)];

    line(H,[0 0],'LineWidth',3,'color','k');
    IQ = quantile(sort(BF),[.25 .75]);
    line(IQ,[0 0],'LineWidth',12,'color','k');
    line([log(3) log(3)],[-.5 .5],'LineWidth',3,'color','r');
    title(['p(logBF>log(3)) = ',num2str(mean(BF>(log(3))))]);
    set(gca,'YTick',[]);
    xlabel('logBF');
    set(gca,'YTick',[]);
    problBFgtlog3(c)=mean(BF>log(3));
    medianlBF(c)=median(BF);
    hdilBF(c,:)=H;
end

out_f=strcat('vbRSA_', outroot,'_', betaroot);

print(out_f,'-dpng','-r300');
save(out_f,'medianlBF','problBFgtlog3','Fs','Fc','hdilBF','b_dist','C','voxels_kept');
%close all;
