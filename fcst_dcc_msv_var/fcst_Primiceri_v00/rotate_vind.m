function new_vind = rotate_vind(vind, mind, nv)
% for "nv" variable models


% rotate_vind so that it matches with mind == 1
% vind is index based on mind == 1
% new_vind is location of the variable based on mind


if nargin == 3 % when nv is specified
    ord_set = flipud(perms([1:nv]));
    cur_set = ord_set(mind,:);
    new_vind = find(cur_set == vind);
    
else %otherwise we assume nv = 3
    
    switch mind
        case 1
            % original ordering
            new_vind = find([1,2,3]==vind);
            
        case 2
            
            new_vind = find([1,3,2]==vind);
        case 3
            
            new_vind = find([2,1,3]==vind);
        case 4
            
            new_vind = find([2,3,1]==vind);
        case 5
            
            new_vind = find([3,1,2]==vind);
        case 6
            
            new_vind = find([3,2,1]==vind);
    end
end