function [rows, cols] = getTerminalSize()
% works in either desktop mode or in terminal
    usingTerminal = ~usejava('desktop');

    % use sensible defaults
    rows = 24;
    cols = 80;

    if (ismac || isunix) && usingTerminal
        % actual terminal: get terminal width using tput
        cmd = 'tput lines';
        [~, r] = unix(cmd);
        num = sscanf(r, '%d');
        if ~isempty(num)
            rows = num;
        end

        cmd = 'tput cols';
        [~, r] = unix(cmd);
        num = sscanf(r, '%d');
        if ~isempty(num)
            cols = num;
        end

    elseif ~usingTerminal
        % matlab command window size
        try
            jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
            cmdWin = jDesktop.getClient('Command Window');
            
            jTextArea = cmdWin.getComponent(0).getViewport.getComponent(0);
            height = jTextArea.getHeight();
            width = jTextArea.getParent.getWidth();
            font = jTextArea.getParent().getFont();
            metrics = cmdWin.getFontMetrics(font);
            charWidth = metrics.charWidth('A');
            charHeight = metrics.getHeight();
            
            rows = floor(height/charHeight);
            cols = floor(width/charWidth);
        catch
        end
    end
end
