classdef (Abstract) FilterClass
   properties
        system
   end
   methods (Abstract)
      prediction(obj,x,y,dt)
   end
end

