function out1 = plotroc(varargin)
%PLOTROC Plot receiver operating characteristic.
%
% <a href="matlab:doc plotroc">plotroc</a>(targets,outputs) takes target data in 1-of-N form (each column
% vector is all zeros with a single 1 indicating the class number), and
% output data and generates a receiver operating characteristic plot.
%
% The best classifications will show the receiver operating line hugging
% the left and top sides of the plots axis.
%
% <a href="matlab:doc plotroc">plotroc</a>(targets,1,outputs1,'name1',targets2,outputs2,names2,...)
% generates a variable number of confusion plots in one figure.
%
% Here a pattern recognition network is trained and its accuracy plotted:
%
%   [x,t] = <a href="matlab:doc simpleclass_dataset">simpleclass_dataset</a>;
%   net = <a href="matlab:doc patternnet">patternnet</a>(10);
%   net = <a href="matlab:doc train">train</a>(net,x,t);
%   y = net(x);
%   <a href="matlab:doc plotroc">plotroc</a>(t,y);
%
% See also roc, plotconfusion, ploterrhist, plotregression.

% Copyright 2010-2014 The MathWorks, Inc.

%% =======================================================
%  BOILERPLATE_START
%  This code is the same for all Transfer Functions.

  persistent INFO;
  if isempty(INFO), INFO = get_info; end
  if nargin == 0
    fig = nnplots.find_training_plot(mfilename);
    if nargout > 0
      out1 = fig;
    elseif ~isempty(fig)
      figure(fig);
    end
    return;
  end
  in1 = varargin{1};
  if ischar(in1)
    switch in1
      case 'info',
        out1 = INFO;
      case 'data_suitable'
        data = varargin{2};
        out1 = nnet.train.isNotParallelData(data);
      case 'suitable'
        [args,param] = nnparam.extract_param(varargin,INFO.defaultParam);
        [net,tr,signals] = deal(args{2:end});
        update_args = standard_args(net,tr,signals);
        unsuitable = unsuitable_to_plot(param,update_args{:});
        if nargout > 0
          out1 = unsuitable;
        elseif ~isempty(unsuitable)
          for i=1:length(unsuitable)
            disp(unsuitable{i});
          end
        end
      case 'training_suitable'
        [net,tr,signals,param] = deal(varargin{2:end});
        update_args = training_args(net,tr,signals,param);
        unsuitable = unsuitable_to_plot(param,update_args{:});
        if nargout > 0
          out1 = unsuitable;
        elseif ~isempty(unsuitable)
          for i=1:length(unsuitable)
            disp(unsuitable{i});
          end
        end
      case 'training'
        [net,tr,signals,param] = deal(varargin{2:end});
        update_args = training_args(net,tr,signals);
        fig = nnplots.find_training_plot(mfilename);
        if isempty(fig)
          fig = figure('Visible','off','Tag',['TRAINING_' upper(mfilename)]);
          plotData = setup_figure(fig,INFO,true);
        else
          plotData = get(fig,'UserData');
        end
        set_busy(fig);
        unsuitable = unsuitable_to_plot(param,update_args{:});
        if isempty(unsuitable)
          set(0,'CurrentFigure',fig);
          plotData = update_plot(param,fig,plotData,update_args{:});
          update_training_title(fig,INFO,tr)
          nnplots.enable_plot(plotData);
        else
          nnplots.disable_plot(plotData,unsuitable);
        end
        fig = unset_busy(fig,plotData);
        if nargout > 0, out1 = fig; end
      case 'close_request'
        fig = nnplots.find_training_plot(mfilename);
        if ~isempty(fig),close_request(fig); end
      case 'check_param'
        out1 = ''; % TODO
      otherwise,
        try
          out1 = eval(['INFO.' in1]);
        catch me, nnerr.throw(['Unrecognized first argument: ''' in1 ''''])
        end
    end
  else
    [args,param] = nnparam.extract_param(varargin,INFO.defaultParam);
    update_args = standard_args(args{:});
    if ischar(update_args)
      nnerr.throw(update_args);
    end
    [plotData,fig] = setup_figure([],INFO,false);
    unsuitable = unsuitable_to_plot(param,update_args{:});
    if isempty(unsuitable)
      plotData = update_plot(param,fig,plotData,update_args{:});
      nnplots.enable_plot(plotData);
    else
      nnplots.disable_plot(plotData,unsuitable);
    end
    set(fig,'Visible','on');
    drawnow;
    if nargout > 0, out1 = fig; end
  end
end

function set_busy(fig)
  set(fig,'UserData','BUSY');
end

function close_request(fig)
  ud = get(fig,'UserData');
  if ischar(ud)
    set(fig,'UserData','CLOSE');
  else
    delete(fig);
  end
  drawnow;
end

function fig = unset_busy(fig,plotData)
  ud = get(fig,'UserData');
  if ischar(ud) && strcmp(ud,'CLOSE')
    delete(fig);
    fig = [];
  else
    set(fig,'UserData',plotData);
  end
  drawnow;
end

function tag = new_tag
  tagnum = 1;
  while true
    tag = [upper(mfilename) num2str(tagnum)];
    fig = nnplots.find_plot(tag);
    if isempty(fig), return; end
    tagnum = tagnum+1;
  end
end

function [plotData,fig] = setup_figure(fig,info,isTraining)
  PTFS = nnplots.title_font_size;
  if isempty(fig)
    fig = get(0,'CurrentFigure');
    if isempty(fig) || strcmp(get(fig,'NextPlot'),'new')
      if isTraining
        tag = ['TRAINING_' upper(mfilename)];
      else
        tag = new_tag;
      end
      fig = figure('Visible','off','Tag',tag);
      if isTraining
        set(fig,'CloseRequestFcn',[mfilename '(''close_request'')']);
      end
    else
      clf(fig);
      set(fig,'Tag','');
      set(fig,'Tag',new_tag);
    end
  end
  set(0,'CurrentFigure',fig);
  ws = warning('off','MATLAB:Figure:SetPosition');
  plotData = setup_plot(fig);
  warning(ws);
  if isTraining
    set(fig,'NextPlot','new');
    update_training_title(fig,info,[]);
  else
    set(fig,'NextPlot','replace');
    set(fig,'Name',[info.name ' (' mfilename ')']);
  end
  set(fig,'NumberTitle','off','ToolBar','none');
  plotData.CONTROL.text = uicontrol('Parent',fig,'Style','text',...
    'Units','normalized','Position',[0 0 1 1],'FontSize',PTFS,...
    'FontWeight','bold','ForegroundColor',[0.7 0 0]);
  set(fig,'UserData',plotData);
end

function update_training_title(fig,info,tr)
  if isempty(tr)
    epochs = '0';
    stop = '';
  else
    epochs = num2str(tr.num_epochs);
    if isempty(tr.stop)
      stop = '';
    else
      stop = [', ' tr.stop];
    end
  end
  set(fig,'Name',['Neural Network Training ' ...
    info.name ' (' mfilename '), Epoch ' epochs stop]);
end

%  BOILERPLATE_END
%% =======================================================

function info = get_info
  info = nnfcnPlot(mfilename,'Receiver Operating Characteristic',7.0,[]);
end

function args = training_args(net,tr,data)
  yall  = nncalc.y(net,data.X,data.Xi,data.Ai);
  y = {yall};
  t = {gmultiply(data.train.mask,data.T)};
  names = {'Training'};
  if ~isempty(data.val.enabled)
    y = [y {yall}];
    t = [t {gmultiply(data.val.mask,data.T)}];
    names = [names {'Validation'}];
  end
  if ~isempty(data.test.enabled)
    y = [y {yall}];
    t = [t {gmultiply(data.test.mask,data.T)}];
    names = [names {'Test'}];
  end
  if length(t) >= 2
    t = [t {data.T}];
    y = [y {yall}];
    names = [names {'All'}];
  end
  args = {t y names};
end

function args = standard_args(varargin)
  if nargin < 2
    args = 'Not enough input arguments.';
  elseif (nargin > 2) && (rem(nargin,3) ~= 0)
    args = 'Incorrect number of input arguments.';
  elseif nargin == 2
    % (t,y)
    t = { nntype.data('format',varargin{1}) };
    y = { nntype.data('format',varargin{2}) };
    names = {''};
    args = {t y names};
  else
    % (t1,y1,name1,...)
    % TODO - Check data is consistent
    count = nargin/3;
    t = cell(1,count);
    y = cell(1,count);
    names = cell(1,count);
    for i=1:count
      t{i} = nntype.data('format',varargin{i*3-2});
      y{i} = nntype.data('format',varargin{i*3-1});
      names{i} = varargin{i*3};
    end
    param.outputIndex = 1;
    args = {t y names};
  end
end

function plotData = setup_plot(fig)
  plotData.numSignals = 0;
end

function fail = unsuitable_to_plot(param,t,y,names)
  fail = '';
  t1 = t{1};
  if numsamples(t1) == 0
    fail = 'The target data has no samples to plot.';
  elseif numtimesteps(t1) == 0
    fail = 'The target data has no timesteps to plot.';
  elseif sum(numelements(t1)) == 0
    fail = 'The target data has no elements to plot.';
  end
end

function plotData = update_plot(param,fig,plotData,tt,yy,names)

  t = tt{1};
  numSignals = length(names);
  numClasses = size(t{1},1);

  % Rebuild figure
  if (plotData.numSignals ~= numSignals) || (plotData.numClasses ~= numClasses)
    set(fig,'NextPlot','replace');
    plotData.numSignals = numSignals;
    plotData.numClasses = numClasses;
    plotData.axes = zeros(1,numSignals);
    colors = nncolor.ncolors(numClasses);
    plotcols = ceil(sqrt(numSignals));
    plotrows = ceil(numSignals/plotcols);
    for plotrow=1:plotrows
      for plotcol=1:plotcols
        i = (plotrow-1)*plotcols+plotcol;
        if (i<=numSignals)
          a = subplot(plotrows,plotcols,i);
          cla(a)
          set(a,'DataAspectRatio',[1 1 1]);
          set(a,'XLim',[0 1]);
          set(a,'YLim',[0 1]);
          hold on
          axisdata = [];
          axisdata.lines = zeros(1,numClasses);
          for j=1:numClasses
            c = colors(j,:);
            line([0 1],[0 1],'LineWidth',2,'Color',[1 1 1]*0.8);
            axisdata.lines(j) = line([0 1],[0 1],'LineWidth',2,'Color',c);
          end
          if ~isempty(names{1})
            titleStr = [names{i} ' ROC'];
          else
            titleStr = 'ROC';
          end
          title(a,titleStr);
          xlabel(a,'False Positive Rate');
          ylabel(a,'True Positive Rate');
          plotData.axes(i) = a;
          set(a,'UserData',axisdata);
          if (i==1) && (numClasses > 1)
            labels = cell(1,numClasses);
            for ii=1:numClasses, labels{ii} = ['Class ' num2str(ii)]; end
            legend(axisdata.lines,labels{:})
          end
        end
      end
    end
    if(~strcmp(fig.WindowStyle, 'docked'))
        screenSize = get(0,'ScreenSize');
        screenSize = screenSize(3:4);
        windowSize = 700 * [1 (plotrows/plotcols)];
        pos = [(screenSize-windowSize)/2 windowSize];
        set(fig,'Position',pos);
    end
  end

  % Update details
  for i=1:numSignals
    y = yy{i}; if iscell(y), y = cell2mat(y); end
    t = tt{i}; if iscell(t), t = cell2mat(t); end
    [tpr,fpr] = roc(t,y);
    if ~iscell(tpr)
      tpr = {tpr};
      fpr = {fpr};
    end
    a = plotData.axes(i);
    axisdata = get(a,'UserData');
    for j=1:numClasses
      set(axisdata.lines(j),'XData',[0 fpr{j} 1],'YData',[0 tpr{j} 1]);
    end
  end
end


