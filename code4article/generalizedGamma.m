%% TODO
% write unittest

classdef generalizedGamma
   properties
      power
      scale
      shift
      %{
      ´p:=power´, ´s:=scale´, ´k:=shift´ are the parameters of the
      (integrand) function:
      ´exp(-x/s)*(x+k)^p´
      %}
   end

   methods

      function obj = generalizedGamma(p,s,k)
         if nargin > 0
            obj.power = p;
            obj.scale = s;
            obj.shift = k;
         end
      end
      
      function estremo_superiore_integrale = evaluate_above(obj,b)
        v = 0;
        for d=0:[obj.power]
            v = v + factorial([obj.power])/factorial(d)*[obj.scale]^([obj.power]-d+1)*( - exp(-b/[obj.scale])*(b+[obj.shift])^(d) );
        end
        estremo_superiore_integrale = v;
      end

      function estremo_inferiore_integrale = evaluate_below(obj,a)
        v = 0;
        for d=0:[obj.power]
            v = v + factorial([obj.power])/factorial(d)*[obj.scale]^([obj.power]-d+1)*( - exp(-a/[obj.scale])*(a+[obj.shift])^(d) );
        end
        estremo_inferiore_integrale = v;
       end

      function integrale = evaluate_integral(obj,a,b)
        integrale = evaluate_above(obj,b) - evaluate_below(obj,a);
      end

   end
end