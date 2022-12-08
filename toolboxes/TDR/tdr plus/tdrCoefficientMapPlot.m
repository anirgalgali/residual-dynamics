function [hf,ha,hi] = tdrCoefficientMapPlot(coefVal,coefChn,coefDes,cNum,cPosx,cPosy,pars)


    

% Number of channels
nch = length(cNum);

% Number of dimensions
ndm = size(coefVal,2);

if isempty(coefDes)
    coefDes = cell(ndm,1);
    for idm = 1:ndm
        coefDes{idm} = sprintf('Dim %i',idm);
    end
end    

% Non-linear transformation
tfun = pars.transform;

% Average within channel
if pars.channelAverage
    cVal = NaN(nch,ndm);
    for ich = 1:nch
        cVal(ich,:) = tfun(mean(coefVal(coefChn==cNum(ich),:)));
    end
    cLoc = cNum;
else
    cVal = tfun(coefVal);
    cLoc = coefChn;
end

% Number of coefficients per location
ncoef = zeros(nch,1);
for ich = 1:nch
    ncoef(ich) = sum(cLoc==cNum(ich));
end

% Set color map
figure; colormap('jet');
jmap = colormap;
cmap = jmap;
cmap(65,:) = [.7 .7 .7];
close gcf;

% Page margins
margX = pars.margin(1);
margY = pars.margin(2);

% Number of segments
nsegX = length(unique(cPosx));
nsegY = length(unique(cPosy));

% Square size at each location
lraw = ceil(sqrt(ncoef));

% Rows and columns
nrow = floor(sqrt(ndm));
ncol = ceil(ndm/nrow);

% Maximum axis size
maxlraw = max(lraw);
maxlnrmX = (1-(ncol+1)*margX)/ncol/nsegX;
maxlnrmY = (1-(nrow+1)*margY)/nrow/nsegY;

% Centers range
crangeX = max(cPosx)-min(cPosx);
crangeY = max(cPosy)-min(cPosy);
minX = min(cPosx);
minY = min(cPosy);

% Normalized position
nrmX = (1-(ncol+1)*margX)/ncol*(cPosx - minX)/crangeX;
nrmY = (1-(nrow+1)*margY)/nrow*(cPosy - minY)/crangeY;

% Find axis for title
candy = find(cPosy==max(cPosy)); 
[~,candx] = min(cPosx(candy));
ides = candy(candx);

% Color limits
for idm = 1:ndm
    switch pars.colorLimit
        case 'full'
            cvlim(idm,:) = [min(cVal(:,idm)) max(cVal(:,idm))];
        case 'abs'
            cvlim(idm,:) = [ 0 +1]*max(abs(cVal(:,idm)));
        case 'sign'
            cvlim(idm,:) = [-1 +1]*max(abs(cVal(:,idm)));
    end
end
mvlim = [min(cvlim(:)) max(cvlim(:))];

% Plot
if isempty(pars.figureHandle)
    hf=figure; axis off;
end
set(hf,'paperposition',[0.25 0.25 8*nrow/ncol 8]);
colormap(cmap);

% Reference dimension
if ~isempty(pars.referenceDimension)    
    
    % The reference
    iref = pars.referenceDimension;
    
    % Reference coefficients
    rr = cVal(:,iref);

    % Loop over channels
    rsort = cell(nch,1);
    for ich = 1:nch
        % Coefficients for this channel
        rrp = rr(cLoc==ich);

        % 1d-sorting
        [~,rsort{ich}] = sort(rrp,'descend');
    end
    
else
    % No 1d-sorting
    rsort = [];
    
end

% Loop over coefficients
hi = zeros(ndm,nch);
ha = zeros(ndm,nch);
for idm = 1:ndm
    
    % Coefficient limits
    switch pars.colorScale
        case 'each'
            vlim = cvlim(idm,:);
        case 'common'
            vlim = mvlim;
    end
    
    % Rows and columns
    icol = mod(idm-1,ncol)+1;
    irow = nrow - floor(idm/ncol-0.5/ncol);
    
    % Coefficients to plot
    vv = cVal(:,idm);
    
    % Loop over channels
    for ich = 1:nch
        
        % Coefficients for this position
        vvp = vv(cLoc==ich);
        
        % The axes center
        axctrX = margX*icol + nrmX(ich) + (icol-1)*(1-margX*(ncol+1))/ncol;
        axctrY = margY*irow + nrmY(ich) + (irow-1)*(1-margY*(nrow+1))/nrow;
        
        % The axes width
        wx = maxlnrmX * lraw(ich) / maxlraw;
        wy = maxlnrmY * lraw(ich) / maxlraw;
        
        % Axes
        ha(idm,ich)=axes('position',[axctrX-wx/2 axctrY-wy/2 wx wy]);
        set(gca,'xtick',[],'ytick',[]);
        
        % Sorting
        if isempty(rsort)
            isort = 1:ncoef(ich);
        else
            isort = rsort{ich};
        end
        
        % Sorted data in 2d-plot
        [jjx,jjy] = meshgrid(1:lraw(ich),1:lraw(ich));
        jjs = jjx+jjy;
        [~,matsort] = sort(jjs(:),'ascend');
        
        % Matrix to plot
        jj = nan(lraw(ich),lraw(ich));
        jj(matsort(1:ncoef(ich))) = vvp(isort);
        
        % Transform into color indeces
        ii = ceil((jj-vlim(1))*(64-1)/diff(vlim)+1);
        ii(isnan(ii)) = 65;
        
        % Plot coefficients
        hi(idm,ich)=image(ii);
        axis off;
        
        % Title
        if ich==ides
            axlim = get(gca,'xlim');
            aylim = get(gca,'ylim');
            ht=text(axlim(1),aylim(1),sprintf('%s    [%.1f  %.1f]',coefDes{idm},vlim(1),vlim(2)));
            set(ht,'horizontalalignment','left','verticalalignment','bottom',...
                'interpreter','none');
        end
        
    end
end


return


%%

% Pool coefficients
cvalue = cell(nax,1);
for iax = 1:nax
    cvalue{iax} = [coef(:,iax).value]';
end
coefVal = [cvalue{:}];

% Pool dsp numbers
coefChn = [coef(:,1).dspnum]';

% Labels
coefDes = [];

% Positions
cPosx = round(1e4 * arrayDspPosition(iar).posX)/1e4;
cPosy = round(1e4 * arrayDspPosition(iar).posY)/1e4;

% Dsp channels
cNum = arrayDspPosition(iar).dspnum;

% Parameters
pars = [];

% Average all coefficients for a given channel
pars.channelAverage = 0;
% pars.channelAverage = 1;

% Nonlinear transformation of coefficient values
% pars.transform = @(x)x;
pars.transform = @(x)sign(x).*log(1+abs(x*50));
% pars.transform = @(x)log(1+abs(x*50));

% Color scaling
% pars.colorScale = 'each';
pars.colorScale = 'common';

% Color limit
% pars.colorLimit = 'full';
% pars.colorLimit = 'abs';
pars.colorLimit = 'sign';

% Ordering within channel
% pars.referenceDimension = [];
pars.referenceDimension = 1;

% Margin
pars.margin = [0.07 0.07];

% Figure
pars.figureHandle = [];
% pars.figureHandle = hf;

% Plot
[hf,ha,hi] = tdrCoefficientMapPlot(coefVal,coefChn,coefDes,cNum,cPosx,cPosy,pars);


