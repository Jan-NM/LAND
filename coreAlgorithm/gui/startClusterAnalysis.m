function varargout = startClusterAnalysis(varargin)
% STARTCLUSTERANALYSIS GUI script for ClusterAnalysis evaluation.
%   
%   
% Cremer Group, Institute of Molecular Biology (IMB), Mainz

% Last Modified by GUIDE v2.5 30-Apr-2019 11:27:15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @startClusterAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @startClusterAnalysis_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%% Function to load content into the listbox 'filebox'
function load_filebox(dir_path, handles)

filter_str = get(handles.fn,'String');   
dir_struct= dir([dir_path '/' filter_str]); 
[sorted_names,sorted_index] = sortrows({dir_struct.name}');

handles.file_names = sorted_names;
handles.is_dir = [dir_struct.isdir];
handles.sorted_index = sorted_index;
guidata(handles.figure1,handles)
set(handles.filebox,'String',handles.file_names,...
 'Value',1)


%% Function to load content into the listbox 'selbox' from the 'filebox'
function load_selbox(fnames, handles)

[~, x] = size(fnames);
[sorted_names,sorted_index] = sortrows(fnames);
handles.file_names = sorted_names;
handles.sorted_index = sorted_index;
guidata(handles.figure1,handles)
set(handles.filenum,'String',[num2str(x) ' files selected'])
set(handles.selbox,'String',handles.file_names,...
 'Value',1)


%% Callback: openingFcn
% --- Executes just before startClusterAnalysis is made visible.
function startClusterAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to startClusterAnalysis (see VARARGIN)

% Choose default command line output for startClusterAnalysis
global files;

DIR_input = pwd;
DIR_output = [pwd filesep 'out'];

files ={};
set(handles.displayDIR_input,'String',DIR_input);
set(handles.displayDIR_output,'String',DIR_output);
load_filebox(DIR_input,handles);

handles.output = hObject;

guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = startClusterAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in selDIR_input.
function selDIR_input_Callback(hObject, eventdata, handles)
DIR_input = get(handles.displayDIR_input,'String');
if(isempty(DIR_input))
    startdir = pwd;
else
    startdir = DIR_input;
end
newDIR_input = uigetdir(startdir,'Select Input Directory');
if(newDIR_input == 0) % user pressed cancel
    return;
else
    DIR_input = newDIR_input;
end
DIR_output = [DIR_input filesep 'out'];
set(handles.filebox,'String',DIR_input);
guidata(hObject, handles);
set(handles.displayDIR_input,'String',DIR_input);
set(handles.displayDIR_output,'String',DIR_output);
load_filebox(DIR_input,handles);


% --- Executes on selection change in filebox.
function filebox_Callback(hObject, eventdata, handles)
global files;
contents = cellstr(get(hObject,'String'));
if(isempty(contents))
    return;
end
newfilename = contents{get(hObject,'Value')};
myfiles = get(handles.selbox,'String');
x = strmatch(newfilename, myfiles);  %% check if file is already selected
if(isempty(x))                  
    files(end+1) ={newfilename};
    load_selbox(files, handles)
end


% --- Executes during object creation, after setting all properties.
function filebox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in selbox.
function selbox_Callback(hObject, eventdata, handles)
global files;
[a b ]= size(files);

contents = cellstr(get(hObject,'String'));
if(isempty(contents))
    return;
end
filenames = contents{get(hObject,'Value')};
newfiles = {};
j=1;
for i = 1:b    
    if (strcmp(files(i),filenames)==0)
        newfiles(j) = files(i);     
        j=j+1;
    end
end
files = newfiles;
load_selbox(files, handles)


% --- Executes during object creation, after setting all properties.
function selbox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function fn_Callback(hObject, eventdata, handles)
DIR_input = get(handles.displayDIR_input,'String');
load_filebox(DIR_input,handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function fn_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in refresh.
function refresh_Callback(hObject, eventdata, handles)
DIR_input = get(handles.displayDIR_input,'String');
load_filebox(DIR_input,handles);


% --- Executes on button press in loadall.
function loadall_Callback(hObject, eventdata, handles)
global files;
files = cellstr(get(handles.filebox,'String'))';
load_selbox(files, handles);


% --- Executes on button press in selDIR_output.
function selDIR_output_Callback(hObject, eventdata, handles)
% global DIR_output
DIR_output = get(handles.displayDIR_output,'String');
newDIR_output = uigetdir('','Select Output Directory');
if(newDIR_output == 0)
    return;
else
    DIR_output = newDIR_output;
end
set(handles.displayDIR_output,'String',DIR_output);
% hObject    handle to selDIR_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% Main routine: START BUTTON PRESSED
% #######################################################
% ###         start main evaluation routine          ####
% #######################################################
% --- Executes on button press in startbut.
function startbut_Callback(hObject, eventdata, handles)

global files

set(handles.startbut,'Enable','off');
drawnow;
set(handles.errors,'String','no Problems so far...');
%--------------------------------------------------------------------------
% initialize callbacks
Clusterparamstruct.displayStatus = @(text) displayStatus(handles.displaystatus,text);
Clusterparamstruct.displayError  = @(text) displayStatus(handles.errors,text);
Clusterparamstruct.enableStartButton = @() enableStartButton(handles.startbut);

% init general parameters
if get(handles.twoDim, 'Value') == 1
    Clusterparamstruct.dimension = 2;
elseif get(handles.threeDim, 'Value') == 1
    Clusterparamstruct.dimension = 3;
end
Clusterparamstruct.compareRandomData = get(handles.compareRandom, 'Value');
if Clusterparamstruct.compareRandomData == 1
    Clusterparamstruct.randomAlgorithm = get(handles.CSR, 'Value');
    if get(handles.avgDensity, 'Value') == 1
        Clusterparamstruct.densityAlgorithm = 'averageDensity';
        Clusterparamstruct.randomAlgorithmValue = []; % should stey empty - parameter is not needed
    end
    if get(handles.kernelDensity, 'Value') == 1
        Clusterparamstruct.densityAlgorithm = 'kernelDensity';
        Clusterparamstruct.randomAlgorithmValue(1) = str2double(get(handles.sampDist, 'String'));
        Clusterparamstruct.randomAlgorithmValue(2) = str2double(get(handles.gridSpacing, 'String'));
    end
    if get(handles.specValue, 'Value') == 1
        Clusterparamstruct.densityAlgorithm = 'specificValue';
        Clusterparamstruct.randomAlgorithmValue = str2double(get(handles.randomValue, 'String'));
    end
else
    Clusterparamstruct.randomAlgorithm = 1; % create random data for object construction even if not requested
    Clusterparamstruct.densityAlgorithm = 'averageDensity';
    Clusterparamstruct.randomAlgorithmValue = [];
end
Clusterparamstruct.showPlots = get(handles.showPlots, 'Value');

% parameter for kNN distance algortihm
Clusterparamstruct.kNNDistance.algorithm = get(handles.kNNDistance, 'Value');
Clusterparamstruct.kNNDistance.k = str2double(get(handles.kValue, 'String'));
Clusterparamstruct.kNNDistance.maxDistance = str2double(get(handles.maxNNDist, 'String')); % in nm

% parameter for DBSCAN algorithm
Clusterparamstruct.DBSCAN.algorithm = get(handles.DBSCANalgo, 'Value');
Clusterparamstruct.DBSCAN.radius = str2double(get(handles.DBSCANradius, 'String'));
Clusterparamstruct.DBSCAN.minPoint = str2double(get(handles.DBSCANpoints, 'String'));
Clusterparamstruct.DBSCAN.maxDiameter = str2double(get(handles.DBSCANdelete, 'String'));

% parameter for distance analysis
Clusterparamstruct.Distance.algorithm = get(handles.distanceAlgo, 'Value');
Clusterparamstruct.Distance.maxDistance = str2double(get(handles.maxDist, 'String')); % in nm

% parameter for RDF
Clusterparamstruct.RDF.algorithm = get(handles.RDF, 'Value');
Clusterparamstruct.RDF.binSize = str2double(get(handles.binSize, 'String')); % in nm
Clusterparamstruct.RDF.maxDistance = str2double(get(handles.maxDistRDF, 'String')); % in nm

% parameter for Ripley's function
Clusterparamstruct.ripley.algorithm = get(handles.ripley, 'Value');
Clusterparamstruct.ripley.radius = str2double(get(handles.ripleyRadius, 'String')); % in nm
Clusterparamstruct.ripley.maxDistance = str2double(get(handles.ripleyDistance, 'String')); % in nm

% get name of sample
if isempty(get(handles.treatmentName, 'String'))
    Clusterparamstruct.treatmentName = 'clusterData';
else
    Clusterparamstruct.treatmentName = get(handles.treatmentName, 'String');
end

set(handles.displaystatus,'String','Evaluation running...');
%--------------------------------------------------------------------------
%% - create outputfolder if not existing
Clusterparamstruct.DIR_input = get(handles.displayDIR_input,'String');
Clusterparamstruct.DIR_output = get(handles.displayDIR_output,'String');
files = cellstr(get(handles.selbox,'String'))';
if(strcmp(get(handles.filenum,'String'),'no files selected'))
    files = {};
end
if(isempty(files))
    set(handles.displaystatus,'String','Error - Evaluation not successful');
    set(handles.errors,'String','No File Selected');
    set(handles.startbut,'Enable','on');
    error('No files selected');
end

if(~isdir(Clusterparamstruct.DIR_output))  
    mkdir(Clusterparamstruct.DIR_output); 
    disp('output folder created')
end
%% - start evaluation - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Clusterparamstruct.handles = handles;
Clusterparamstruct.outputmode = 'gui';

clusterBatchProcessing(Clusterparamstruct);

set(handles.displaystatus,'String','Evaluation finished!');
set(handles.startbut,'Enable','on');


% --- Executes when selected object is changed in filetypegroup.
function filetypegroup_SelectionChangeFcn(hObject, eventdata, handles)
switch get(eventdata.NewValue,'Tag') 
    case 'ortesel'
        set(handles.fn,'String','*.mat');     
end
DIR_input = get(handles.displayDIR_input,'String');
load_filebox(DIR_input,handles);
% hObject    handle to the selected object in filetypegroup 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function filetypegroup_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to filetypegroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function filetypegroup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filetypegroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function uimenuEnableStartButton_Callback(hObject, eventdata, handles)
% hObject    handle to uimenuEnableStartButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.startbut,'Enable','on');


% --- Executes on button press in debugbutton.
function debugbutton_Callback(hObject, eventdata, handles)
% hObject    handle to debugbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard

 
function displayStatus(hObject,text)
%    call this function with hObject = handles.displaystatus
    set(hObject,'String',text);


function enableStartButton(hObject)
%    call this function with hObject = handles.startbut
    set(hObject,'Enable','on');


function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox10.
function listbox10_Callback(hObject, eventdata, handles)
% hObject    handle to listbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox10


% --- Executes during object creation, after setting all properties.
function listbox10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function order_Callback(hObject, eventdata, handles)
% hObject    handle to order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of order as text
%        str2double(get(hObject,'String')) returns contents of order as a double


% --- Executes during object creation, after setting all properties.
function order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in cropOrte.
function cropOrte_Callback(hObject, eventdata, handles)
% hObject    handle to cropOrte (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cropImage;


% --- Executes on button press in helpbutton.
function helpbutton_Callback(hObject, eventdata, handles)
% hObject    handle to helpbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
winopen('help/LAND_manual.pdf');


% --- Executes on button press in twoDim.
function twoDim_Callback(hObject, eventdata, handles)
% hObject    handle to twoDim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.threeDim, 'Value', 0);
    
% Hint: get(hObject,'Value') returns toggle state of twoDim


% --- Executes on button press in threeDim.
function threeDim_Callback(hObject, eventdata, handles)
% hObject    handle to threeDim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.twoDim, 'Value', 0);

% Hint: get(hObject,'Value') returns toggle state of threeDim


% --- Executes on button press in compareRandom.
function compareRandom_Callback(hObject, eventdata, handles)
% hObject    handle to compareRandom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(get(hObject,'Value'))
     set(handles.CSR, 'Enable', 'on');   
     set(handles.avgDensity, 'Enable', 'on');
     set(handles.kernelDensity, 'Enable', 'on');
     set(handles.specValue, 'Enable', 'on');
     set(handles.sampDist, 'Enable', 'on');
     set(handles.randomValue, 'Enable', 'on');
     set(handles.gridSpacing, 'Enable', 'on');
     set(handles.text80, 'Enable', 'on');
     set(handles.text81, 'Enable', 'on');
else
    set(handles.CSR, 'Enable', 'off');
    set(handles.avgDensity, 'Enable', 'off');
    set(handles.kernelDensity, 'Enable', 'off');
    set(handles.specValue, 'Enable', 'off');
    set(handles.sampDist, 'Enable', 'off');
    set(handles.randomValue, 'Enable', 'off');
    set(handles.gridSpacing, 'Enable', 'off')
    set(handles.text80, 'Enable', 'off');
    set(handles.text81, 'Enable', 'off');
end
% Hint: get(hObject,'Value') returns toggle state of compareRandom


% --- Executes on button press in showPlots.
function showPlots_Callback(hObject, eventdata, handles)
% hObject    handle to showPlots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showPlots



function kValue_Callback(hObject, eventdata, handles)
% hObject    handle to kValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kValue as text
%        str2double(get(hObject,'String')) returns contents of kValue as a double


% --- Executes during object creation, after setting all properties.
function kValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in kNNDistance.
function kNNDistance_Callback(hObject, eventdata, handles)
% hObject    handle to kNNDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of kNNDistance



function maxNNDist_Callback(hObject, eventdata, handles)
% hObject    handle to maxNNDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxNNDist as text
%        str2double(get(hObject,'String')) returns contents of maxNNDist as a double


% --- Executes during object creation, after setting all properties.
function maxNNDist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxNNDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DBSCANradius_Callback(hObject, eventdata, handles)
% hObject    handle to DBSCANradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DBSCANradius as text
%        str2double(get(hObject,'String')) returns contents of DBSCANradius as a double


% --- Executes during object creation, after setting all properties.
function DBSCANradius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DBSCANradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DBSCANpoints_Callback(hObject, eventdata, handles)
% hObject    handle to DBSCANpoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DBSCANpoints as text
%        str2double(get(hObject,'String')) returns contents of DBSCANpoints as a double


% --- Executes during object creation, after setting all properties.
function DBSCANpoints_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DBSCANpoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DBSCANdelete_Callback(hObject, eventdata, handles)
% hObject    handle to DBSCANdelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DBSCANdelete as text
%        str2double(get(hObject,'String')) returns contents of DBSCANdelete as a double


% --- Executes during object creation, after setting all properties.
function DBSCANdelete_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DBSCANdelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DBSCANalgo.
function DBSCANalgo_Callback(hObject, eventdata, handles)
% hObject    handle to DBSCANalgo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DBSCANalgo


% --- Executes on button press in distanceAlgo.
function distanceAlgo_Callback(hObject, eventdata, handles)
% hObject    handle to distanceAlgo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of distanceAlgo



function maxDist_Callback(hObject, eventdata, handles)
% hObject    handle to maxDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxDist as text
%        str2double(get(hObject,'String')) returns contents of maxDist as a double


% --- Executes during object creation, after setting all properties.
function maxDist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in parameterEstimation.
function parameterEstimation_Callback(hObject, eventdata, handles)
% hObject    handle to parameterEstimation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.mat', 'Select an Orte file to test DBSCAN parameter');
if isequal(filename, 0)
    disp('No Orte file selected')
else
    estimationFile = load([pathname filename]);
    oldField = char(fieldnames(estimationFile));
    if ~(strcmp(oldField, 'Orte'))
        newField = 'Orte';
        [estimationFile.(newField)] = estimationFile.(oldField);
        estimationFile = rmfield(estimationFile, oldField);   
    end
    parameterObject = ClusterAnalysis(estimationFile.Orte);
    parameterObject.parameterEstimation('DBSCAN');
    parameterObject.delete;
    clear('Orte', 'parameterObject', 'filename', 'pathname');
end


% --- Executes on button press in avgDensity.
function avgDensity_Callback(hObject, eventdata, handles)
% hObject    handle to avgDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.kernelDensity, 'Value', 0);
set(handles.specValue, 'Value', 0);
% Hint: get(hObject,'Value') returns toggle state of avgDensity


% --- Executes on button press in kernelDensity.
function kernelDensity_Callback(hObject, eventdata, handles)
% hObject    handle to kernelDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.avgDensity, 'Value', 0);
set(handles.specValue, 'Value', 0);
% Hint: get(hObject,'Value') returns toggle state of kernelDensity


% --- Executes on button press in specValue.
function specValue_Callback(hObject, eventdata, handles)
% hObject    handle to specValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.avgDensity, 'Value', 0);
set(handles.kernelDensity, 'Value', 0);
% Hint: get(hObject,'Value') returns toggle state of specValue


function sampDist_Callback(hObject, eventdata, handles)
% hObject    handle to sampDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sampDist as text
%        str2double(get(hObject,'String')) returns contents of sampDist as a double


% --- Executes during object creation, after setting all properties.
function sampDist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sampDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function randomValue_Callback(hObject, eventdata, handles)
% hObject    handle to randomValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of randomValue as text
%        str2double(get(hObject,'String')) returns contents of randomValue as a double


% --- Executes during object creation, after setting all properties.
function randomValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to randomValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CSR.
function CSR_Callback(hObject, eventdata, handles)
% hObject    handle to CSR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CSR



function treatmentName_Callback(hObject, eventdata, handles)
% hObject    handle to treatmentName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of treatmentName as text
%        str2double(get(hObject,'String')) returns contents of treatmentName as a double


% --- Executes during object creation, after setting all properties.
function treatmentName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to treatmentName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function binSize_Callback(hObject, eventdata, handles)
% hObject    handle to binSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of binSize as text
%        str2double(get(hObject,'String')) returns contents of binSize as a double


% --- Executes during object creation, after setting all properties.
function binSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to binSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RDF.
function RDF_Callback(hObject, eventdata, handles)
% hObject    handle to RDF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RDF



function maxDistRDF_Callback(hObject, eventdata, handles)
% hObject    handle to maxDistRDF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxDistRDF as text
%        str2double(get(hObject,'String')) returns contents of maxDistRDF as a double


% --- Executes during object creation, after setting all properties.
function maxDistRDF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxDistRDF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ripleyRadius_Callback(hObject, eventdata, handles)
% hObject    handle to ripleyRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ripleyRadius as text
%        str2double(get(hObject,'String')) returns contents of ripleyRadius as a double


% --- Executes during object creation, after setting all properties.
function ripleyRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ripleyRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ripleyDistance_Callback(hObject, eventdata, handles)
% hObject    handle to ripleyDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ripleyDistance as text
%        str2double(get(hObject,'String')) returns contents of ripleyDistance as a double


% --- Executes during object creation, after setting all properties.
function ripleyDistance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ripleyDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ripley.
function ripley_Callback(hObject, eventdata, handles)
% hObject    handle to ripley (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ripley



function gridSpacing_Callback(hObject, eventdata, handles)
% hObject    handle to gridSpacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gridSpacing as text
%        str2double(get(hObject,'String')) returns contents of gridSpacing as a double


% --- Executes during object creation, after setting all properties.
function gridSpacing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gridSpacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
