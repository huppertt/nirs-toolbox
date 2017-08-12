function disable_plot(plotData,unsuitable)

% Copyright 2010 The MathWorks, Inc.
  string = [{'','',''},unsuitable];
  set(plotData.CONTROL.text,'string',string,'visible','on');
end
