function [ XX ] = cut_data( X, percent_top, percent_bottom, dir )
% Removes data from top/bottom along axis
% dir should be a row vector

    if (percent_top > 1)
        percent_top = percent_top/100;
    end
    if (percent_bottom > 1)
        percent_bottom = percent_bottom/100;
    end
    
    Y = X*(dir');
    minY = min(Y);
    maxY = max(Y);
    len = maxY-minY;
    
    c_bottom = len*percent_bottom + minY;
    c_top = maxY - len*percent_top;
    
    XX = X((Y >= c_bottom) & (Y <= c_top), :);

end

