function h=veloViewer(array,array2,voxelDim,cmap)
% arrayviewer - arrayviewer for looking at 3d or 4d date in gui
% array = magnitude
% array2 - phase
%
%      usage: [  ] = veloViewer(array,array2,tr,cmap)
%         by: jonas larsson
%        mod: alex beckett
%       date: 2015-07-02
%        $Id: veloViewer.m 
%     inputs: array, array2, voxelDim, cmap
%    outputs: h (figure handle)
%
%    purpose: program originally written by jonas larsson to have a quick look at
%    3d/4d data without having to load things into complicated
%    programs.
%
%    modified by alex beckett to quickly observe mag/phase data from PC imaging
%
%        e.g: arrayviewer(rand(20,30,40,50), jet)
%

% initialize and check data.
hdr=[];
if (isstruct(array))
    % TODO -- need to add some logic here to check!
    hdr=array.hdr;
    array=array.data;
end
if ~isnumeric(array) && ~islogical(array)
    disp('(arrayviewer) Array must be numeric')
    return;
end
if ~isfloat(array)
    array = single(array);
end



nd=ndims(array);
%if (nd<3 || nd>4)
sizeArray = size(array);
if nd == 2
    sizeArray(3:4) = 1;
elseif nd == 3
    sizeArray(4) = 1;
end
if nd>4
    error('array must be 3d or 4d')
end

array=imrotate(array,-90);
array2=imrotate(array2,-90);

data.array=mean(array,4);
data.array2=array2;

% figure out limits to show
data = getDataLimits(data);
if isempty(data.limits)
    return
end

% otherwise make figure
h=figure;


data.hdr=hdr;
data.currorient=[1 2 3];
data.currslicedim=[3 1 2];
data.roiselect=0;
data.currslice=round(sizeArray/2);
if (nargin>=3) && length(voxelDim)==4;
    data.voxelDim=voxelDim;
    
else
    data.voxelDim=[1 1 1 1];
end
data.tr=data.voxelDim(4);
if (nargin==4)
    data.cmap=cmap;
else
    data.cmap=gray;
end
set(h,'KeyPressFcn',@keypress);
%set(h,'WindowButtonDownFcn',@buttonpress);
set(h,'WindowScrollWheelFcn',@scroll);
set(h,'UserData',data);

drawXXslice(h);

end

function keypress(src,evnt)
currdata=get(src,'UserData');
currdata.roiselect=0;
switch (evnt.Key)
    case {'h','H','a','A','1'}
        currdata.currorient=1;
        currdata.currslicedim=3;
    case {'s','S','2'}
        currdata.currorient=2;
        currdata.currslicedim=1;
    case {'c','C','3'}
        currdata.currorient=3;
        currdata.currslicedim=2;
    case 'pageup'
        if (length(currdata.currorient)>1)
            return
        end
        currdata.currslice(currdata.currslicedim)=...
            currdata.currslice(currdata.currslicedim)+10;
    case 'pagedown'
        if (length(currdata.currorient)>1)
            return
        end
        currdata.currslice(currdata.currslicedim)=...
            currdata.currslice(currdata.currslicedim)-10;
    case 'uparrow'
        if (length(currdata.currorient)>1)
            return
        end
        currdata.currslice(currdata.currslicedim)=...
            currdata.currslice(currdata.currslicedim)+1;
    case 'downarrow'
        if (length(currdata.currorient)>1)
            return
        end
        currdata.currslice(currdata.currslicedim)=...
            currdata.currslice(currdata.currslicedim)-1;
    case 'home'
        if (length(currdata.currorient)>1)
            return
        end
        currdata.currslice(currdata.currslicedim)=round(size(currdata.array,currdata.currslicedim)/2);
    case 'leftarrow'
        currdata.currslice(4)=currdata.currslice(4)-1;
    case 'rightarrow'
        currdata.currslice(4)=currdata.currslice(4)+1;
    case {'o','4'}
        currdata.currorient=[1 2 3];
        currdata.currslicedim=[3 1 2];
        currdata.roiselect=0;
    case {'x','5'}
        currdata.currorient=[1 2 3];
        currdata.currslicedim=[3 1 2];
        currdata.roiselect=1;
    otherwise
        return
end
set(src,'UserData',currdata);
drawXXslice(src);
end

function scroll(src,evnt)
currdata=get(src,'UserData');
if (length(currdata.currorient)>1)
    return
end
currdata.currslice(currdata.currslicedim)=...
    currdata.currslice(currdata.currslicedim)+evnt.VerticalScrollCount;
set(src,'UserData',currdata);
drawXXslice(src);

end
function buttonpress(src,evnt)
prnt=get(src,'Parent');
prntprnt=get(prnt,'Parent');
currdata=get(prntprnt,'UserData');
p=round(get(prnt,'CurrentPoint'));
whichorient=get(prnt,'UserData');
p=p(1,1:2);
if (any(p<1))
    return
end

P=[1 1 1]';
switch (whichorient)
    case 1
        if (p(1)>size(currdata.array,1) || p(2)>size(currdata.array,2))
            return
        end
        P=[p(1,1),p(1,2),currdata.currslice(3)]';
    case 2
        if (p(1)>size(currdata.array,2) || p(2)>size(currdata.array,3))
            return
        end
        P=[currdata.currslice(1),p(1,1),p(1,2)]';
    case 3
        if (p(1)>size(currdata.array,1) || p(2)>size(currdata.array,3))
            return
        end
        P=[p(1,1),currdata.currslice(2),p(1,2)]';
    case 4
        if (p(1)>size(currdata.array,4))
            return
        end
        P=currdata.currslice(1:3)';
        currdata.currslice(4)=p(1);
end
currdata.currslice(1:3)=P;
set(prntprnt,'UserData',currdata);
drawXXslice(prntprnt)

end

function drawslice(h)
% drawslice function --- call back
figure(h)
currdata=get(h,'UserData');
plotorder=[2 1 4];
if (currdata.currslice(4)<1)
    currdata.currslice(4)=1;
elseif (currdata.currslice(4)>size(currdata.array,4))
    currdata.currslice(4)=size(currdata.array,4);
end
for n=1:length(currdata.currorient)
    if (currdata.currslice(currdata.currslicedim(n))<1)
        currdata.currslice(currdata.currslicedim(n))=1;
    elseif (currdata.currslice(currdata.currslicedim(n))>...
            size(currdata.array,currdata.currslicedim(n)))
        currdata.currslice(currdata.currslicedim(n))=size(currdata.array,currdata.currslicedim(n));
    end
    if (length(currdata.currorient)>1)
        currdata.subplots(n)=subplot(2,2,plotorder(n));
    else
        if (size(currdata.array2,4)>1)
            currdata.subplots=subplot(2,1,1);
        else
            currdata.subplots=subplot(1,1,1);
        end
    end
    switch (currdata.currorient(n))
        case 1
            img = currdata.array(:,:,currdata.currslice(3),currdata.currslice(4))';
            x=[.5 size(currdata.array,1)+.5; currdata.currslice(1)  currdata.currslice(1)]';
            y=[currdata.currslice(2)  currdata.currslice(2); .5 size(currdata.array,2)+.5]';
        case 2
            img = squeeze(currdata.array(currdata.currslice(1),:,:,currdata.currslice(4)))';
            x=[.5 size(currdata.array,2)+.5; currdata.currslice(2)  currdata.currslice(2)]';
            y=[currdata.currslice(3)  currdata.currslice(3); .5 size(currdata.array,3)+.5]';
        case 3
            img = squeeze(currdata.array(:,currdata.currslice(2),:,currdata.currslice(4)))';
            x=[.5 size(currdata.array,1)+.5; currdata.currslice(1)  currdata.currslice(1)+.5]';
            y=[currdata.currslice(3)  currdata.currslice(3); .5 size(currdata.array,3)+.5]';
    end
    
    isNanOrInf = isnan(img) | isinf(img);
    if any(isNanOrInf(:))
        [nanXCoords,nanYCoords] = ind2sub(size(img),find(isNanOrInf));
        img(isNanOrInf) = size(currdata.cmap,1);
        imagesc(img,currdata.limits');
        hold on
        numberNanVoxels = length(nanXCoords);
        %draw a contour around NaNbox centered on zero
        nanVoxelX = [-.5 .5;-.5 .5;-.5 -.5;.5 .5];%;-.5 .5];
        nanVoxelY = [-.5 -.5;.5 .5;-.5 .5;-.5 .5];%;.5 -.5];
        nanVoxelCoords = [nanVoxelX nanVoxelY];
        %replicate this pattern at all nan positions
        nanXCoords = reshape(repmat(nanXCoords',numel(nanVoxelX),1),2,size(nanVoxelX,1)*numberNanVoxels)';
        nanYCoords = reshape(repmat(nanYCoords',numel(nanVoxelY),1),2,size(nanVoxelY,1)*numberNanVoxels)';
        nanCoords = [nanXCoords nanYCoords]+repmat(nanVoxelCoords,numberNanVoxels,1);
        %remove common segments (segments that appear twice)
        [~, uniqueSegments] = unique(nanCoords,'rows');
        commonSegments = setdiff(1:size(nanCoords,1),uniqueSegments);
        %nanCoords(commonSegments,:)=[];
        nanCoords = setdiff(nanCoords,nanCoords(commonSegments,:),'rows');
        uniqueSegmentCoords = unique(nanCoords,'rows');
        commonSegmentCoords = setdiff(nanCoords,uniqueSegmentCoords,'rows');
        nanCoords = setdiff(nanCoords,commonSegmentCoords,'rows');
        nanXCoords = nanCoords(:,1:2)';
        nanYCoords = nanCoords(:,3:4)';
        plot(nanYCoords,nanXCoords,'k'); %swith X and Y
    else
        imagesc(img,  currdata.limits);%liyong chen 2013-06-12
        hold on
    end
    
    colormap(currdata.cmap)
    %colorbar
    %axis image
    axis xy
    plot(x,y,'r')
    drawnow
    
    hold off
    set(get(currdata.subplots(n),'Children'),'ButtonDownFcn',@buttonpress);
    set(currdata.subplots(n),'UserData',currdata.currorient(n));
end
if ( size(currdata.array2,4)>1)
    if (length(currdata.currorient)>1)
        ch=subplot(2,2,3);
    else
        ch=subplot(2,1,2);
    end
    plot(squeeze(currdata.array2(currdata.currslice(1),currdata.currslice(2),currdata.currslice(3),:)));
    hold on
    plot([currdata.currslice(4) currdata.currslice(4)],get(ch,'ylim'),'r')
    set(gca, 'xtick', unique( [get(ch,'xlim') currdata.currslice(4) ]))
    hold off
    set(get(ch,'Children'),'ButtonDownFcn',@buttonpress);
    set(ch,'UserData',4);
    
end
set(h,'UserData',currdata);

%set the figure name
P = currdata.currslice;
v=currdata.array(P(1),P(2),P(3),P(4));
if (isempty(currdata.hdr) || ~isfield(currdata.hdr,'qform44'))
    s=sprintf('Voxel at %i %i %i, value=%f\n',P(1:3),v);
else
    QP=currdata.hdr.qform44*[P;1];
    SP=currdata.hdr.sform44*[P;1];
    s=sprintf('Voxel at %i %i %i\nScanner coords %03.2f %03.2f %03.2f\nAligned coords %03.2f %03.2f %03.2f\nValue=%f\n',P,QP(1:3),SP(1:3),v);
end
set(h,'Name',s)

end

function [data] = getDataLimits(data)
% helper function to get data limits (and avoid nan's inf's in pathological
% cases,,,

% in particular remove nans
% and infs.
robustRange = [5 95];
array = data.array(~isinf(data.array) & ~isnan(data.array));
if isempty(array)
    disp('(getDataLimits) There are only nans/infs in this array')
    data.limits = [];
    return
end

data.limits= prctile(array(:), robustRange); % [min(array(:)) max(array(:))];

if data.limits(1)==data.limits(2)%if that returns two identical values
    %use nanmin and nanmax
    data.limits(1) = nanmin(array(:));
    data.limits(2) = nanmax(array(:));
end

if data.limits(1)==data.limits(2)% || ~issorted(data.limits) %can't remember what that's for
    fprintf('(getDataLimits) there is only one value in the array (%d)',data.limits(1));
    if ~data.limits(1)
        data.limits = [-1 1];
    else
        data.limits = data.limits+.1*[-1 1].*abs(data.limits);
    end
end

end




function drawXXslice(h)
% drawslice function --- call back
figure(h)
currdata=get(h,'UserData');
plotorder=[2 1 4];
if (currdata.currslice(4)<1)
    currdata.currslice(4)=1;
elseif (currdata.currslice(4)>size(currdata.array,4))
    currdata.currslice(4)=size(currdata.array,4);
end
for n=1:length(currdata.currorient)
    if (currdata.currslice(currdata.currslicedim(n))<1)
        currdata.currslice(currdata.currslicedim(n))=1;
    elseif (currdata.currslice(currdata.currslicedim(n))>...
            size(currdata.array,currdata.currslicedim(n)))
        currdata.currslice(currdata.currslicedim(n))=size(currdata.array,currdata.currslicedim(n));
    end
    if (length(currdata.currorient)>1)
        currdata.subplots(n)=subplot(2,2,plotorder(n));
    else
        if (size(currdata.array2,4)>1)
            currdata.subplots=subplot(2,1,1);
        else
            currdata.subplots=subplot(1,1,1);
        end
    end
    switch (currdata.currorient(n))
        case 1
            img = currdata.array(:,:,currdata.currslice(3),currdata.currslice(4))';
            x=[.5 size(currdata.array,1)+.5; currdata.currslice(1)  currdata.currslice(1)]';
            y=[currdata.currslice(2)  currdata.currslice(2); .5 size(currdata.array,2)+.5]';
        case 2
            img = squeeze(currdata.array(currdata.currslice(1),:,:,currdata.currslice(4)))';
            x=[.5 size(currdata.array,2)+.5; currdata.currslice(2)  currdata.currslice(2)]';
            y=[currdata.currslice(3)  currdata.currslice(3); .5 size(currdata.array,3)+.5]';
        case 3
            img = squeeze(currdata.array(:,currdata.currslice(2),:,currdata.currslice(4)))';
            x=[.5 size(currdata.array,1)+.5; currdata.currslice(1)  currdata.currslice(1)+.5]';
            y=[currdata.currslice(3)  currdata.currslice(3); .5 size(currdata.array,3)+.5]';
    end
    
    isNanOrInf = isnan(img) | isinf(img); 
    if any(isNanOrInf(:))
        [nanXCoords,nanYCoords] = ind2sub(size(img),find(isNanOrInf));
        img(isNanOrInf) = size(currdata.cmap,1);
        imagesc(img,currdata.limits);
        hold on
        numberNanVoxels = length(nanXCoords);
        %draw a contour around NaNbox centered on zero
        nanVoxelX = [-.5 .5;-.5 .5;-.5 -.5;.5 .5];%;-.5 .5];
        nanVoxelY = [-.5 -.5;.5 .5;-.5 .5;-.5 .5];%;.5 -.5];
        nanVoxelCoords = [nanVoxelX nanVoxelY];
        %replicate this pattern at all nan positions
        nanXCoords = reshape(repmat(nanXCoords',numel(nanVoxelX),1),2,size(nanVoxelX,1)*numberNanVoxels)';
        nanYCoords = reshape(repmat(nanYCoords',numel(nanVoxelY),1),2,size(nanVoxelY,1)*numberNanVoxels)';
        nanCoords = [nanXCoords nanYCoords]+repmat(nanVoxelCoords,numberNanVoxels,1);
        %remove common segments (segments that appear twice)
        [~, uniqueSegments] = unique(nanCoords,'rows');
        commonSegments = setdiff(1:size(nanCoords,1),uniqueSegments);
        %nanCoords(commonSegments,:)=[];
        nanCoords = setdiff(nanCoords,nanCoords(commonSegments,:),'rows');
        uniqueSegmentCoords = unique(nanCoords,'rows');
        commonSegmentCoords = setdiff(nanCoords,uniqueSegmentCoords,'rows');
        nanCoords = setdiff(nanCoords,commonSegmentCoords,'rows');
        nanXCoords = nanCoords(:,1:2)';
        nanYCoords = nanCoords(:,3:4)';
        plot(nanYCoords,nanXCoords,'k'); %swith X and Y
    else
        imagesc(img,  currdata.limits);%liyong chen 2013-06-12
        daspect(1./currdata.voxelDim([setdiff(1:3,currdata.currslicedim(n)),currdata.currslicedim(n)]));
        hold on
    end
    
    colormap(currdata.cmap)
    %colorbar
    %axis image
    axis xy
    plot(x,y,'r')
   
    drawnow
    hold off
    set(get(currdata.subplots(n),'Children'),'ButtonDownFcn',@buttonpress);
    set(currdata.subplots(n),'UserData',currdata.currorient(n));
end
if ( size(currdata.array2,4)>1)
    
    %timelength=size(currdata.array2,4);
    if(currdata.roiselect<1)
        if (length(currdata.currorient)>1)
            ch=subplot(2,2,3);
        else
            ch=subplot(2,1,2);
        end
        plotData=(squeeze(currdata.array2(currdata.currslice(1),currdata.currslice(2),currdata.currslice(3),:)));
        Fs=1/currdata.tr;
        d=fdesign.bandpass('N,F3dB1,F3dB2',10,.8,1.4,Fs);
        %Hd=design(d,'butter');
        %filtData=filtfilthd(Hd,plotData);
        d=fdesign.lowpass('Fp,Fst,Ap,Ast',.2,.5,2,60,Fs);
        %Hd=design(d,'butter');
        %filtData_Resp=filtfilthd(Hd,plotData);
        plot(plotData,'k-')
        line([0 length(plotData)],[0 0],'linestyle','--','color','k');

        hold on
        %plot(filtData,'r')
        %plot(filtData_Resp,'b')
        %set(gca,'xTick',[0:30:(timelength)])

        %set(gca,'xTickLabel',[(0:30:(timelength))*currdata.tr])
%         set(gca,'FontSize',14)
         %xlabel('Time(unit:s)')%,'fontsize',22)
         %ylabel('Velocity(unit:cm/s)')%,'fontsize',22)
         %title('CSF velocity')%,'fontsize',22)
%         hold off
        axis('tight');
        %ylim([-pi pi]);
        xlim([1 length(plotData)]);
    elseif 0
        figure(20)
        imgtmpm=currdata.array(:,:,currdata.currslice(3));
        maxPtV=max(reshape(imgtmpm,1,[]));
        BW=roipoly(imgtmpm/(maxPtV/2));
        %BW_all{pivot_at}=BW;
        for kk=1:size(currdata.array2,4)
            imgtmp=squeeze(currdata.array2(:,:,currdata.currslice(3),kk));
            imgp(kk)=mean(imgtmp(BW));
        end
        figure(h)
        if (length(currdata.currorient)>1)
            ch=subplot(2,2,3);
        else
            ch=subplot(2,1,2);
        end
        
        plot(imgp*2.5/pi)
        set(gca,'xTick',[0:50:(timelength-1)])
        set(gca,'xTickLabel',[(0:50:(timelength-1))*.0825])
        set(gca,'FontSize',14)
        xlabel('Time(unit:s)','fontsize',22)
        ylabel('Velocity(unit:cm/s)','fontsize',22)
        title('CSF velocity','fontsize',22)
        
        subplot(2,2,2)
        %imagesc
        imgtmpm(BW)=1000;
        imagesc(imgtmpm',currdata.limits)
        hold off
%         hold on
%         x=[.5 size(currdata.array,1)+.5; currdata.currslice(1)  currdata.currslice(1)]';
%         y=[currdata.currslice(2)  currdata.currslice(2); .5 size(currdata.array,2)+.5]';
%         axis image
%         axis xy
%         plot(x,y,'r')
        %plot(squeeze(currdata.array2(currdata.currslice(1),currdata.currslice(2),currdata.currslice(3),:)));
    end
    hold on
    %plot([currdata.currslice(4) currdata.currslice(4)],get(ch,'ylim'),'r')
    %set(gca, 'xtick', unique( [get(ch,'xlim') currdata.currslice(4) ]))
    hold off
    set(get(ch,'Children'),'ButtonDownFcn',@buttonpress);
    set(ch,'UserData',4);
    if 0
        figure(20)
        subplot(1,2,1)
        plot(plotData,'k')
        set(gca,'xTick',[0:30:(timelength)])
            set(gca,'xTickLabel',[(0:30:(timelength))*currdata.tr])
            set(gca,'FontSize',14)
            xlabel('Time(unit:s)')
            ylabel('Velocity(unit:cm/s)')
            title('CSF velocity')
            axis('tight')





    %     subplot(2,2,3)
    %     plot(filtData)
    %     set(gca,'xTick',[0:50:(timelength-1)])
    %         set(gca,'xTickLabel',[(0:50:(timelength-1))*currdata.tr])
    %         set(gca,'FontSize',14)
    %         xlabel('Time(unit:s)')
    %         ylabel('Velocity(unit:cm/s)')
    %         title('CSF velocity, Band-Pass Filtered (Cardiac)')
        Fs=1/currdata.tr;
        T=currdata.tr;
        L=length(plotData);
        t=(0:L-1)*T;
        NFFT=2^nextpow2(L);
        Y=fft(plotData,NFFT)/L;
        f=Fs/2*linspace(0,1,NFFT/2+1);
        subplot(1,2,2)
        plot(f,2*abs(Y(1:NFFT/2+1)));
        axis('tight');
            title('Single-Sided Amplitude Spectrum of y(t)');xlabel('Frequency(Hz)');ylabel('|Y(f)|')
        try 
            load('physioVar')

            resp=respSR{currdata.currslice(3),1};
            for i=1:size(resp,1)
                B(i,:)=pinv([resp(i,:)' ones(size(resp,2),1)])*plotData;
            end
            [~,j]=max(B(:,1));
            subplot(1,2,1)
            hold on
            plot([resp(j,:)' ones(size(resp,2),1)]*B(j,:)','b--')
            subplot(1,2,2)
            hold on
            ft=fft([resp(j,:)' ones(size(resp,2),1)]*B(j,:)',NFFT)/L;
            ft(1)=0;
            plot(f,abs(ft(1:NFFT/2+1)),'b--');
        catch
            disp('No physio data in current directory')
        end
    end
%     subplot(2,2,4)
%     plot(filtData_Resp)
%     set(gca,'xTick',[0:50:(timelength-1)])
%         set(gca,'xTickLabel',[(0:50:(timelength-1))*currdata.tr])
%         set(gca,'FontSize',14)
%         xlabel('Time(unit:s)')
%         ylabel('Velocity(unit:cm/s)')
%         title('CSF velocity, Lo-Pass Filtered (Resp)')
    
    
end
set(h,'UserData',currdata);

%set the figure name
P = currdata.currslice;
v=currdata.array(P(1),P(2),P(3),P(4));
if (isempty(currdata.hdr) || ~isfield(currdata.hdr,'qform44'))
    s=sprintf('Voxel at %i %i %i, value=%f\n',P(1:3),v);
else
    QP=currdata.hdr.qform44*[P;1];
    SP=currdata.hdr.sform44*[P;1];
    s=sprintf('Voxel at %i %i %i\nScanner coords %03.2f %03.2f %03.2f\nAligned coords %03.2f %03.2f %03.2f\nValue=%f\n',P,QP(1:3),SP(1:3),v);
end
set(h,'Name',s)

end


function x=filtfilthd(varargin)
% FILTFILTHD Zero-phase digital filtering with dfilt objects.
%
% FILTFILLTHD provides zero phase filtering and accepts dfilt objects on
% input. A number of end-effect minimization methods are supported.
%
% Examples:
% x=FILTFILTHD(Hd, x)
% x=FILTFILTHD(Hd, x, method)
% where Hd is a dfilt object and x is the input data. If x is a matrix,
% each column will be filtered.
%
% ------------------------------------------------------------------------
% The filter states of Hd on entry to FILTFILTHD will be used at the
% beginning of each forward and each backward pass through the data. The 
% user should normally ensure that the initial states are zeroed before
% calling filtfilthd [e.g. using reset(Hd);]
% ------------------------------------------------------------------------
%
% x=FILTFILTHD(b, a, x) 
% x=FILTFILTHD(b, a, x, method)
%           format is also supported.
%
% x=FILTFILTHD(...., IMPULSELENGTH)
%   allows the impulse response length to be specified on input. 
%
%
% method is a string describing the end-effect correction technique:
%   reflect:      data at the ends are reflected and mirrored as in
%                   the MATLAB filtfilt function (default)
%   predict:      data are extraploated using linear prediction
%   spline/pchip: data are extrapolated using MATLAB's interp1 function
%   none:         no internal end-effect correction is applied
%                   (x may be pre-pended and appended with data externally)
%
% Each method has different merits/limitations. The most robust
% method is reflect.
% 
% The length of the padded data section at each end will be impzlength(Hd)
% points, or with 'reflect', the minimum of impzlength(Hd) and the
% data length (this is different to filtfilt where the padding is only
% 3 * the filter width). Using the longer padding reduces the need for any
% DC correction (see the filtfilt documentation).
%
%
% See also: dfilt, filtfilt, impzlength
%
% -------------------------------------------------------------------------
% Author: Malcolm Lidierth 10/07
% Copyright � The Author & King's College London 2007
% -------------------------------------------------------------------------
%
% Revisions:
%   07.11.07    nfact=len-1 (not len) when impulse response longer than
%                           data
%   11.11.07    Use x=f(x) for improved memory performance 
%                           [instead of y=f(x)]
%               Handles row vectors properly
%   31.01.08    Allow impzlength to be specified on input
%                           

% ARGUMENT CHECKS
if isnumeric(varargin{1}) && isnumeric(varargin{2})
    % [b, a] coefficients on input, convert to a singleton filter.
    Hd = dfilt.df2(varargin{1}, varargin{2});
    x=varargin{3};
elseif ~isempty(strfind(class(varargin{1}),'dfilt'))
    % dfilt object on input
    Hd=varargin{1};
    x=varargin{2};
else
    error('Input not recognized');
end

if ischar(varargin{end})
    method=varargin{end};
else
    method='reflect';
end

if isscalar(varargin{end})
    nfact=varargin{end};
else
    nfact=impzlength(Hd);
end

% DEAL WITH MATRIX INPUTS-----------------------------------------------
% Filter each column in turn through recursive calls
[m,n]=size(x);
if (m>1) && (n>1) 
    for i=1:n
        x(:,i)=filtfilthd(Hd, x(:,i), method, nfact);
    end
    return
end

% MAIN FUNCTION-------------------------------------------------------
% Make sure x is a column. Return to row vector later if needed
if m==1
    x = x(:);
    trflag=true;
else
    trflag=false;
end

len=length(x);
switch method
    case 'reflect'
        % This is similar to the MATLAB filtfilt reflect and mirror method
        nfact=min(len-1, nfact);%change to len-1 not len 07.11.07
        pre=2*x(1)-x(nfact+1:-1:2);
        post=2*x(len)-x(len-1:-1:len-nfact);
    case 'predict'
        % Use linear prediction. DC correction with mean(x).
        % Fit over 2*nfact points, with nfact coefficients
        np=2*nfact;
        m=mean(x(1:np));
        pre=lpredict(x(1:np)-m, nfact, nfact, 'pre')+m;
        m=mean(x(end-np+1:end));
        post=lpredict(x(end-np+1:end)-m, nfact, nfact, 'post')+m;
    case {'spline', 'pchip'}
        % Spline/pchip extrapolation.
        % Fit over 2*nfact points,
        np=2*nfact;
        pre=interp1(1:np, x(np:-1:1), np+1:np+nfact, method, 'extrap');
        pre=pre(end:-1:1)';
        post=interp1(1:np, x(end-np+1:end), np+1:np+nfact, method, 'extrap')';
    case 'none'
        % No end-effect correction
        pre=[];
        post=[];
end

% % UNCOMMENT TO VIEW THE PADDED DATA 
% if length(pre);plot(pre,'color', 'r');end;
% line(length(pre)+1:length(pre)+length(x), x);
% if length(post);line(length(pre)+length(x)+1:length(pre)+length(x)+length(post), post, 'color', 'm');end;


% Remember Hd is passed by reference - save entry state and restore later
memflag=get(Hd, 'persistentmemory');
states=get(Hd, 'States');
set(Hd,'persistentmemory', true);

%--------------------------------
% ----- FORWARD FILTER PASS -----
% User-supplied filter states at entry will be applied
% Pre-pended data
pre=filter(Hd, pre); %#ok<NASGU>
% Input data
x=filter(Hd ,x);
% Post-pended data
post=filter(Hd, post);

% ------ REVERSE FILTER PASS -----
% Restore user-supplied filter states for backward pass
set(Hd, 'States', states);
% Post-pended reversed data
post=filter(Hd, post(end:-1:1)); %#ok<NASGU>
% Reversed data
x=filter(Hd, x(end:-1:1));

% Restore data sequence
x=x(end:-1:1);
%---------------------------------
%---------------------------------

% Restore Hd
set(Hd, 'States', states);
set(Hd,'persistentmemory', memflag)

% Revert to row if necessary
if trflag
    x=x.';   
end

return
end


%--------------------------------------------------------------------------
function y=lpredict(x, np, npred, pos)
% LPREDICT estimates the values of a data set before/after the observed
% set.
% LPREDICT Local version. For a stand-alone version see:
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=16798&objectType=FILE
% -------------------------------------------------------------------------
% Author: Malcolm Lidierth 10/07
% Copyright � The Author & King's College London 2007
% -------------------------------------------------------------------------
% Order input sequence
if nargin==4 && strcmpi(pos,'pre')
    x=x(end:-1:1);
end
% Get the forward linear predictor coefficients via the LPC
% function
a=lpc(x,np);
% Negate coefficients, and get rid of a(1)
cc=-a(2:end);
% Pre-allocate output
y=zeros(npred,1);
% Seed y with the first value
y(1)=cc*x(end:-1:end-np+1);
% Next np-1 values
for k=2:min(np,npred)
    y(k)=cc*[y(k-1:-1:1); x(end:-1:end-np+k)];
end
% Now do the rest
for k=np+1:npred
    y(k)=cc*y(k-1:-1:k-np);
end
% Order the output sequence if required
if nargin==4 && strcmpi(pos,'pre')
    y=y(end:-1:1);
end
return
end
% -------------------------------------------------------------------------
