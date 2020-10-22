function h=arrayviewer(array,cmap,limits)
% arrayviewer -  h=arrayviewer(array,cmap,limits)
%
%      usage: [  ] = arrayviewer( array, cmap  )
%         by: jonas larsson
%       date: 2009-07-02
%        $Id: arrayviewer.m 830 2011-05-31 07:42:16Z lpzjb $:
%     inputs: array, cmap
%    outputs: h (figure handle)
%
%    purpose: program written by jonas larsson to have a quick look at
%    3d/4d data without having to load things into complicated
%    programs.
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
data.array=array;

if (nargin==3)
    data.limits=limits;
else
    % figure out limits to show
    data = getDataLimits(data);
    if isempty(data.limits)
        return
    end
end

% otherwise make figure
h=figure;

data.hdr=hdr;
data.currorient=1;
data.currslicedim=3;
data.currorient=[1 2 3];
data.currslicedim=[3 1 2];
data.currslice=round(sizeArray/2);
if (nargin>=2)
    data.cmap=cmap;
else
    data.cmap=gray;
end
set(h,'KeyPressFcn',@keypress);
%set(h,'WindowButtonDownFcn',@buttonpress);
set(h,'WindowScrollWheelFcn',@scroll);
set(h,'UserData',data);

drawslice(h)

end

function keypress(src,evnt)
currdata=get(src,'UserData');
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
    otherwise
        return
end
set(src,'UserData',currdata);
drawslice(src);
end

function scroll(src,evnt)
currdata=get(src,'UserData');
if (length(currdata.currorient)>1)
    return
end
currdata.currslice(currdata.currslicedim)=...
    currdata.currslice(currdata.currslicedim)+evnt.VerticalScrollCount;
set(src,'UserData',currdata);
drawslice(src);

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
drawslice(prntprnt)

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
        if (size(currdata.array,4)>1)
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
      [dump, uniqueSegments] = unique(nanCoords,'rows');
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
      imagesc(img,  currdata.limits);
      hold on
    end
    
    colormap(currdata.cmap)
    colorbar
    axis image
    axis xy
    plot(x,y,'r')
    drawnow
    hold off
    set(get(currdata.subplots(n),'Children'),'ButtonDownFcn',@buttonpress);
    set(currdata.subplots(n),'UserData',currdata.currorient(n));
end
if ( size(currdata.array,4)>1)
    if (length(currdata.currorient)>1)
        ch=subplot(2,2,3);
    else
        ch=subplot(2,1,2);
    end
    plot(squeeze(currdata.array(currdata.currslice(1),currdata.currslice(2),currdata.currslice(3),:)));
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
    disp(sprintf('(getDataLimits) there is only one value in the array (%d)',data.limits(1)));
    if ~data.limits(1)
      data.limits = [-1 1];
    else
      data.limits = data.limits+.1*[-1 1].*abs(data.limits);
    end
  end

end

