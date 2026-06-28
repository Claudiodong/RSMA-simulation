function y = squeeze_size(x)
         % check the input matrix size
         [varargout] = size(x);
         N = length(varargout);

         for i = 1:N
               if (varargout(i) == 1)
                   y = squeeze(x);
               end

         end

end