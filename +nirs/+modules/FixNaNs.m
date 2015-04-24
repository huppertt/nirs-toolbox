classdef FixNaNs < nirs.modules.AbstractModule

    properties
        ifFailReplaceWith = 0;
    end
    
    methods

        function obj = FixNaNs( prevJob )
           obj.name = 'Fix NaNs';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:length(data)

                d = data(i).data;
                t = data(i).time;
                
                lst = isnan(d);
                
                if any(lst(:))
                
                    try
                        for j = 1:size(d,2)
                            if any(lst(:,j))
                                % interpolation
                                l = lst(:,j);
                                d(l,j) = interp1(t(~l), d(~l,j), t(l),'linear','extrap');
                                data(i).data = d;
                            end
                        end
                        
                    catch
                        % just replace with zeros
                        d(lst) = obj.ifFailReplaceWith;
                    end

                end
            end
        end
    end
    
end

