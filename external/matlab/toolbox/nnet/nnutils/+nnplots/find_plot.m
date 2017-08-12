function fig = find_plot(tag)

% Copyright 2010 The MathWorks, Inc.
  for object = get(0,'children')'
    if strcmp(get(object,'type'),'figure') 
      if strcmp(get(object,'tag'),tag)
       fig = object;
       return
     end
    end
  end
  fig = [];
end
