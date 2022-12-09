classdef (Abstract) DynamicSys
   properties
      n
      m
      sigmaX
      sigmaY
   end
   methods (Abstract)
      stateTransition(obj,x,dt)
      measurement(obj,x)
   end
end