function setLineOpacity(s, alpha)
% stores information in UserData struct to cause saveFigure to render
% lines as translucent when exporting to svg

error('not supported');
    [edge, face] = cellvec(numel(s));
    for i = 1:numel(s)
        
        % tag it as translucent for saveFigure to
        % pick up during SVG authoring
        userdata = get(s(i),'UserData');
        userdata.svg.LineAlpha = alpha;
        set(s(i),'UserData', userdata);
            
        if ~verLessThan('matlab', '8.4')
            % first cache marker opacity
            if isempty(s(i).MarkerHandle) && ~isa(s(i).MarkerHandle, 'matlab.graphics.GraphicsPlaceholder')
                edge{i} = s(i).MarkerHandle.EdgeColorData;
                face{i} = s(i).MarkerHandle.FaceColorData;
            end
            
            % use RGBA color specification
            s(i).Color(4) = alpha;
        end
    end
    
    drawnow;
    
    for i = 1:length(s)
        if ~isempty(edge{i})
            s(i).MarkerHandle.EdgeColorData = edge;
            s(i).MarkerHandle.FaceColorData = face;
        end
    end

end