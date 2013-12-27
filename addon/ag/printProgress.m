function printProgress(message, progress, maximum, showTime)
    persistent previousLineLength;
    persistent t;
    
    if isempty(previousLineLength) || progress == 0
        previousLineLength = 0;
        t = cputime;
    end
    
    percentage = round(progress * 100 / maximum);
    if exist('showTime', 'var') && showTime
        activity = ['(' num2str(cputime - t) ' s)'];
    else
        animation = '-\|/';
        animationIndex = 1 + mod(progress, length(animation));
        activity = animation(animationIndex);
    end
    
    line = sprintf('%s:% 4.3g %% %s', message, percentage, activity);
    
    fprintf([repmat('\b', [1 previousLineLength]) '%s\n'], line);
    
    previousLineLength = length(line) + 1;
end
