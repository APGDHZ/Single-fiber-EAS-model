function [sortArray,varargout] = sortAll(sortArray,direction,varargin)
    
    if nargin>1
        [sortArray, ind] = sort(sortArray,2,direction) ;
    else
        sortArray = sort(sortArray,2) ;
    end
    for i=1:numel(varargin)
        if min(size(varargin{i})) == 1
            varargout{i} = varargin{i}(ind) ; % sort array
        else
            varargout{i} = varargin{i}(:,ind) ; % sort matrix
        end
    end
    
end

