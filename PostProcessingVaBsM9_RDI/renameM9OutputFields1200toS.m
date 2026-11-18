function [DepthRNGout,WaterProfileOut,BTvelOut,ShearStressOut, ...
    BedloadTransportOut,RenameReport] = renameM9OutputFields1200toS( ...
    DepthRNGout,WaterProfileOut,BTvelOut,ShearStressOut, ...
    BedloadTransportOut,varargin)
%RENAMEM9OUTPUTFIELDS1200TOS Rename RDI 1200/12 field tokens to M9 _S.
%
%   [DepthRNGout,WaterProfileOut,BTvelOut,ShearStressOut, ...
%    BedloadTransportOut,RenameReport] = renameM9OutputFields1200toS( ...
%       DepthRNGout,WaterProfileOut,BTvelOut,ShearStressOut, ...
%       BedloadTransportOut)
%
% replaces both the full instrument token 1200 and the isolated abbreviated
% token 12 with _S in structure field names. An underscore immediately
% before either token is reused, so the result never gets a double
% underscore. The isolated-12 rule does not alter 12 when it is part of a
% larger number such as 0125. For example:
%
%       Depth1200       -> Depth_S
%       Depth_1200      -> Depth_S
%       Vmod_1200       -> Vmod_S
%       RNG11_1200_     -> RNG11_S_
%       BTmodf_1200ALL  -> BTmodf_SALL
%       WF_12           -> WF_S
%       WT12            -> WT_S
%       alpha12mean     -> alpha_Smean
%       Ks2_12          -> Ks2_S
%       QPer0125D       -> QPer0125D  (unchanged)
%
% The function changes names only. Numeric values, strings, filenames,
% units, arrays, and time-series cells are not modified. Nested structures
% (including inherited output inside BedloadTransportOut), structures held
% in cells, and table/timetable variable names are processed recursively.
%
% Name-value options:
%   'Recursive'             true (default) or false
%   'RenameTableVariables'  true (default) or false
%   'CollisionPolicy'       What to do if both a 1200/12 field and its S
%                           target already exist:
%                           'prefer_existing_s' (default)
%                           'prefer_renamed_1200'
%                           'error'
%   'AssignToCaller'        false (default). If true, the five converted
%                           outputs and RenameReport are assigned in the
%                           calling workspace.
%   'Verbose'               true (default) or false
%
% RenameReport records every renamed or removed field/table variable.

p = inputParser;
p.FunctionName = mfilename;
isStructOrEmpty = @(x) isempty(x) || isstruct(x);
addRequired(p,'DepthRNGout',isStructOrEmpty);
addRequired(p,'WaterProfileOut',isStructOrEmpty);
addRequired(p,'BTvelOut',isStructOrEmpty);
addRequired(p,'ShearStressOut',isStructOrEmpty);
addRequired(p,'BedloadTransportOut',isStructOrEmpty);
addParameter(p,'Recursive',true,@(x) islogical(x) && isscalar(x));
addParameter(p,'RenameTableVariables',true,@(x) islogical(x) && isscalar(x));
addParameter(p,'CollisionPolicy','prefer_existing_s', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p,'AssignToCaller',false,@(x) islogical(x) && isscalar(x));
addParameter(p,'Verbose',true,@(x) islogical(x) && isscalar(x));
parse(p,DepthRNGout,WaterProfileOut,BTvelOut,ShearStressOut, ...
    BedloadTransportOut,varargin{:});
opt = p.Results;
opt.CollisionPolicy = lower(string(opt.CollisionPolicy));

validPolicies = ["prefer_existing_s","prefer_renamed_1200","error"];
if ~any(opt.CollisionPolicy == validPolicies)
    error('renameM9OutputFields1200toS:InvalidCollisionPolicy', ...
        ['CollisionPolicy must be ''prefer_existing_s'', ' ...
         '''prefer_renamed_1200'', or ''error''.']);
end

rootNames = ["DepthRNGout","WaterProfileOut","BTvelOut", ...
    "ShearStressOut","BedloadTransportOut"];
values = {DepthRNGout,WaterProfileOut,BTvelOut,ShearStressOut, ...
    BedloadTransportOut};
records = emptyRenameRecords();

for k = 1:numel(values)
    if isempty(values{k}) && ~isstruct(values{k})
        continue
    end
    [values{k},records] = renameValue1200toS(values{k},rootNames(k), ...
        rootNames(k),opt,records);
end

DepthRNGout          = values{1};
WaterProfileOut      = values{2};
BTvelOut             = values{3};
ShearStressOut       = values{4};
BedloadTransportOut  = values{5};

if isempty(records)
    RenameReport = table('Size',[0 6], ...
        'VariableTypes',repmat("string",1,6), ...
        'VariableNames',{'RootStructure','ContainerPath','OldName', ...
        'NewName','ObjectType','Action'});
else
    RenameReport = struct2table(records);
end

remaining = strings(0,1);
for k = 1:numel(values)
    remaining = [remaining; collectLegacyInstrumentNames( ...
        values{k},rootNames(k),opt)]; %#ok<AGROW>
end
remaining = unique(remaining,'stable');
if ~isempty(remaining)
    warning('renameM9OutputFields1200toS:NamesRemain', ...
        ['Some field/table names still contain a 1200 or isolated 12 ', ...
         'instrument token after conversion:\n%s'], ...
        strjoin(remaining,newline));
end

if opt.AssignToCaller
    assignin('caller','DepthRNGout',DepthRNGout);
    assignin('caller','WaterProfileOut',WaterProfileOut);
    assignin('caller','BTvelOut',BTvelOut);
    assignin('caller','ShearStressOut',ShearStressOut);
    assignin('caller','BedloadTransportOut',BedloadTransportOut);
    assignin('caller','RenameReport',RenameReport);
end

if opt.Verbose
    if height(RenameReport) == 0
        fprintf('%s: no names containing 1200 or isolated 12 were found.\n', ...
            mfilename);
    else
        fprintf('%s: processed %d field/table-variable name(s).\n', ...
            mfilename,height(RenameReport));
        nCollision = sum(contains(RenameReport.Action,"collision"));
        if nCollision > 0
            fprintf('  Resolved %d name collision(s) using policy "%s".\n', ...
                nCollision,opt.CollisionPolicy);
        end
    end
end
end


function [value,records] = renameValue1200toS(value,rootName,path,opt,records)
% Dispatch conversion according to the container type.

if isstruct(value)
    [value,records] = renameStruct1200toS(value,rootName,path,opt,records);
elseif istable(value) || isa(value,'timetable')
    if opt.RenameTableVariables
        [value,records] = renameTable1200toS(value,rootName,path,opt,records);
    end
elseif iscell(value) && opt.Recursive
    for n = 1:numel(value)
        itemPath = path + "{" + n + "}";
        [value{n},records] = renameValue1200toS(value{n},rootName, ...
            itemPath,opt,records);
    end
end
end


function [S,records] = renameStruct1200toS(S,rootName,path,opt,records)
% Recursively process structure contents, then rename its fields.

names = fieldnames(S);
if isempty(names)
    return
end

if isempty(S)
    [S,records] = renameEmptyStruct1200toS(S,rootName,path,opt,records);
    return
end

if opt.Recursive
    for f = 1:numel(names)
        thisName = names{f};
        for n = 1:numel(S)
            if numel(S) == 1
                childPath = path + "." + string(thisName);
            else
                childPath = path + "(" + n + ")." + string(thisName);
            end
            [S(n).(thisName),records] = renameValue1200toS( ...
                S(n).(thisName),rootName,childPath,opt,records);
        end
    end
end

names = fieldnames(S);
for f = 1:numel(names)
    oldName = names{f};
    if ~isLegacyInstrumentName(oldName) || ~isfield(S,oldName)
        continue
    end
    newName = convert1200TokenToS(oldName);
    oldPath = path + "." + string(oldName);

    if isfield(S,newName)
        if structFieldsEqual(S,oldName,newName)
            S = rmfield(S,oldName);
            action = "identical_collision_removed_1200";
        else
            switch opt.CollisionPolicy
                case "prefer_existing_s"
                    S = rmfield(S,oldName);
                    action = "collision_kept_existing_s";
                case "prefer_renamed_1200"
                    for n = 1:numel(S)
                        S(n).(newName) = S(n).(oldName);
                    end
                    S = rmfield(S,oldName);
                    action = "collision_replaced_s_with_1200";
                otherwise
                    error('renameM9OutputFields1200toS:FieldCollision', ...
                        ['Cannot rename %s to %s because the target field ' ...
                         'already exists and the values differ.'], ...
                        oldPath,path + "." + string(newName));
            end
        end
    else
        for n = 1:numel(S)
            S(n).(newName) = S(n).(oldName);
        end
        S = rmfield(S,oldName);
        action = "renamed";
    end

    records = appendRenameRecord(records,rootName,path,oldName,newName, ...
        "struct_field",action);
end
end


function [S,records] = renameEmptyStruct1200toS(S,rootName,path,opt,records)
% Rename fields of a structure array with zero elements.

oldNames = fieldnames(S);
finalNames = oldNames;

for f = 1:numel(oldNames)
    oldName = oldNames{f};
    if ~isLegacyInstrumentName(oldName) || ~any(strcmp(finalNames,oldName))
        continue
    end
    newName = convert1200TokenToS(oldName);
    targetExists = any(strcmp(finalNames,newName));

    if targetExists
        if opt.CollisionPolicy == "error"
            error('renameM9OutputFields1200toS:FieldCollision', ...
                ['Cannot rename %s.%s to %s.%s because the target field ' ...
                 'already exists.'],path,oldName,path,newName);
        end
        finalNames(strcmp(finalNames,oldName)) = [];
        if opt.CollisionPolicy == "prefer_renamed_1200"
            action = "empty_collision_preferred_1200_name";
        else
            action = "empty_collision_kept_existing_s";
        end
    else
        finalNames{strcmp(finalNames,oldName)} = newName;
        action = "renamed_empty_struct_field";
    end

    records = appendRenameRecord(records,rootName,path,oldName,newName, ...
        "struct_field",action);
end

if ~isequal(oldNames,finalNames)
    oldSize = size(S);
    template = cell2struct(cell(numel(finalNames),1),finalNames(:),1);
    S = repmat(template,oldSize);
end
end


function [T,records] = renameTable1200toS(T,rootName,path,opt,records)
% Rename table/timetable variables and recurse into their contents.

names = T.Properties.VariableNames;

if opt.Recursive
    for f = 1:numel(names)
        thisName = names{f};
        childPath = path + "." + string(thisName);
        [T.(thisName),records] = renameValue1200toS(T.(thisName), ...
            rootName,childPath,opt,records);
    end
end

names = T.Properties.VariableNames;
for f = 1:numel(names)
    oldName = names{f};
    if ~isLegacyInstrumentName(oldName) || ...
            ~ismember(oldName,T.Properties.VariableNames)
        continue
    end
    newName = convert1200TokenToS(oldName);
    targetExists = ismember(newName,T.Properties.VariableNames);

    if targetExists
        if isequaln(T.(oldName),T.(newName))
            T(:,{oldName}) = [];
            action = "identical_collision_removed_1200";
        else
            switch opt.CollisionPolicy
                case "prefer_existing_s"
                    T(:,{oldName}) = [];
                    action = "collision_kept_existing_s";
                case "prefer_renamed_1200"
                    T.(newName) = T.(oldName);
                    T(:,{oldName}) = [];
                    action = "collision_replaced_s_with_1200";
                otherwise
                    error('renameM9OutputFields1200toS:TableVariableCollision', ...
                        ['Cannot rename table variable %s.%s to %s.%s ' ...
                         'because the target already exists and differs.'], ...
                        path,oldName,path,newName);
            end
        end
    else
        index = strcmp(T.Properties.VariableNames,oldName);
        T.Properties.VariableNames{index} = newName;
        action = "renamed";
    end

    records = appendRenameRecord(records,rootName,path,oldName,newName, ...
        "table_variable",action);
end
end


function tf = structFieldsEqual(S,fieldA,fieldB)
% Compare two fields across every element in a structure array.

tf = true;
for n = 1:numel(S)
    if ~isequaln(S(n).(fieldA),S(n).(fieldB))
        tf = false;
        return
    end
end
end


function newName = convert1200TokenToS(oldName)
% Convert 1200 and isolated 12 tokens to one normalized _S marker.

newName = regexprep(oldName,'_?1200','_S');
newName = regexprep(newName,'_?(?<!\d)12(?!\d)','_S');
end


function tf = isLegacyInstrumentName(name)
% True for 1200 or an isolated 12 token, but false for numbers like 0125.

tf = contains(name,'1200') || ...
    ~isempty(regexp(name,'(?<!\d)12(?!\d)','once'));
end


function records = emptyRenameRecords()
records = struct('RootStructure',{},'ContainerPath',{},'OldName',{}, ...
    'NewName',{},'ObjectType',{},'Action',{});
end


function records = appendRenameRecord(records,rootName,path,oldName, ...
    newName,objectType,action)

entry.RootStructure = string(rootName);
entry.ContainerPath = string(path);
entry.OldName = string(oldName);
entry.NewName = string(newName);
entry.ObjectType = string(objectType);
entry.Action = string(action);
records(end+1) = entry; %#ok<AGROW>
end


function namesFound = collectLegacyInstrumentNames(value,path,opt)
% Find names still containing a 1200 or isolated 12 instrument token.

namesFound = strings(0,1);

if isstruct(value)
    fields = fieldnames(value);
    for f = 1:numel(fields)
        thisName = fields{f};
        if isLegacyInstrumentName(thisName)
            namesFound(end+1,1) = path + "." + string(thisName); %#ok<AGROW>
        end
        if opt.Recursive
            for n = 1:numel(value)
                if numel(value) == 1
                    childPath = path + "." + string(thisName);
                else
                    childPath = path + "(" + n + ")." + string(thisName);
                end
                namesFound = [namesFound; collectLegacyInstrumentNames( ...
                    value(n).(thisName),childPath,opt)]; %#ok<AGROW>
            end
        end
    end
elseif (istable(value) || isa(value,'timetable')) && opt.RenameTableVariables
    variables = value.Properties.VariableNames;
    for f = 1:numel(variables)
        thisName = variables{f};
        if isLegacyInstrumentName(thisName)
            namesFound(end+1,1) = path + "." + string(thisName); %#ok<AGROW>
        end
        if opt.Recursive
            namesFound = [namesFound; collectLegacyInstrumentNames( ...
                value.(thisName), ...
                path + "." + string(thisName),opt)]; %#ok<AGROW>
        end
    end
elseif iscell(value) && opt.Recursive
    for n = 1:numel(value)
        namesFound = [namesFound; collectLegacyInstrumentNames(value{n}, ...
            path + "{" + n + "}",opt)]; %#ok<AGROW>
    end
end
end
