function myCopyObj(axFrom, axTo)
    % Custom version of copyobj() which works to copy a UIAxes object to
    % an old-school axes object.
    
    % Copy children (lines).
    copyobj(axFrom.Children, axTo);
    % Copy titles and labels.
    copyobj([axFrom.Title axFrom.XLabel axFrom.YLabel], axTo)
    
    % Transfer other properties.
    uiAxParams = get(axFrom);
    uiAxParamNames = fieldnames(uiAxParams);
    % Get list of editable params in new axis
    editableParams = fieldnames(set(axTo));
    % Remove the UIAxes params that aren't editable in the new axes (add others you don't want)
    badFields = uiAxParamNames(~ismember(uiAxParamNames, editableParams));
    badFields = [badFields; 'Parent'; 'Children'; 'XAxis'; 'YAxis'; 'ZAxis'; 'XLabel'; 'YLabel'; 'ZLabel'; 'Title'];
    %'Position'; 'OuterPosition'; 'XLabel'; 'YLabel'; 'ZLabel'; 'Title'];
    uiAxGoodParams = rmfield(uiAxParams,badFields);
    % set editable params on new axes
    set(axTo, uiAxGoodParams)
end