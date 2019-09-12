function varargout = AClusterTool(varargin)
% AClusterTool MATLAB code for AClusterTool.fig
%      AClusterTool, by itself, creates a new AClusterTool or raises the existing
%      singleton*.
%
%      H = AClusterTool returns the handle to a new AClusterTool or the handle to
%      the existing singleton*.
%
%      AClusterTool('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AClusterTool.M with the given input arguments.
%
%      AClusterTool('Property','Value',...) creates a new AClusterTool or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AClusterTool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AClusterTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AClusterTool

% Last Modified by GUIDE v2.5 12-Sep-2019 15:57:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AClusterTool_OpeningFcn, ...
                   'gui_OutputFcn',  @AClusterTool_OutputFcn, ...
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


% --- Executes just before AClusterTool is made visible.
function AClusterTool_OpeningFcn(hObject,~, handles, varargin)%eventdata
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AClusterTool (see VARARGIN)

% Choose default command line output for AClusterTool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AClusterTool wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AClusterTool_OutputFcn(~, ~, handles) %hObject, eventdata
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%Aquí se encuentran las funciones para cada botón%%%%

% --- Executes on button press in MatrixSelectionButton.
function MatrixSelectionButton_Callback(~, ~, handles) %hObject, eventdata
% hObject    handle to MatrixSelectionButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[matrix] = uigetfile ('*.mat'); %Standard open file dialog box for .mat files
load([matrix]);
matrix = open(matrix); %Open file in appropriate application.
name = fieldnames(matrix); %Get structure field names or Java/COM object properties.
matrix = matrix.(name{1}); %Set name as the ID of matrix values
matrix = matrix'; %For this event we have to transpose the matrix. Found a way to identify when to do it autom
handles.output = matrix; %Set the matrix obtained as the new handles value
assignin('base','matrix',handles.output);

%Lets have a quick look at how the data look with a histogram
length (Spikes); %Substitue this value in the whole script!!!
poblationV = zeros(1,length(Spikes)); % Preallocation of a matrix 1xN, where N is the length of the matrix (timeFrames)
column = (1); %Set the initial columnID number as 1
t = [1:1:length(Spikes)];

while column <= length(Spikes) %While the actual column doesnt surpass the N of columns existent
    for i=1:length(Spikes) %For every item from 1 to the number of columns
        poblationV(1,column) = sum(Spikes(:,column)); %Select its correspondent place in prealloc matrix
        column = column+1; %Add 1 to the columnID every iteration so we can continue to the next column
        %disp(poblationV);
    end
end

plot (t,poblationV,'.');


% --- Executes on button press in PCAButton.
%Principal Component Analysis
function PCAButton_Callback(hObject, ~, handles)
matrix=evalin('base','matrix');
[X] = (matrix);
no_dims = round(intrinsic_dim(X, 'MLE'));
disp(['MLE estimate of intrinsic dimensionality: ' num2str(no_dims)]);
[mappedX, mapping] = compute_mapping(X, 'PCA', no_dims);
scatter(mappedX(:,1), mappedX(:,2)); title('Result of PCA');
%%matrix=evalin('base','matrix');
%[coeff,score,explained] = pca (matrix);
%scatter3(score(:,1),score(:,2),score(:,3)) %Scatterplot of the 3 first principal components
%axis equal
%xlabel('1st Principal Component')
%ylabel('2nd Principal Component')
%zlabel('3rd Principal Component')
%%guidata(hObject,handles);


% --- Executes on button press in KmeansButton.
function KmeansButton_Callback(hObject, eventdata, handles)
% hObject    handle to KmeansButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mappedX=evalin('base','mappedX');
idx=kmeans(mappedX(:,:),3);
assignin('base','idx',idx);
hold on
    xlabel('Component1'),ylabel('Component2'),zlabel('PCA3')
    colores= [.5 1 .5;1 0 0;0 0 1];%[75 0 130; 32 178 170]./255;
    for n=1:size(mappedX,1)
        plot(mappedX(n,1),mappedX(n,2),'.'...
            ,'color',colores(idx(n),:)) 
    end


    
% --- Executes on button press in MeanShiftButton.
function MeanShiftButton_Callback(hObject, eventdata, handles)
% hObject    handle to MeanShiftButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 matrix=evalin('base','matrix');
    [matrix, mapping] = isomap(matrix, 2, []);
    assignin('base','mappedX',matrix);


% --- Executes on button press in MDSButton.
%Classical Multidimensional Scaling
function MDSButton_Callback(hObject, eventdata, handles)
% hObject    handle to MDSButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
matrix=evalin('base','matrix');
[X] = (matrix);
no_dims = round(intrinsic_dim(X, 'MLE'));
mappedX = mds(X, no_dims);
scatter(mappedX(:,1), mappedX(:,2)); title('Result of MDS');
assignin('base','mappedX',mappedX);
