clear;
clc;

%% About the units
% m day kPa

TIME_STEPS = 3600;
dt = 360;
print_dt = 360;

eventTimes = unique([0:dt:TIME_STEPS, TIME_STEPS]);

%% PATH TO INIT STATE
CurrDir = fileparts(mfilename('fullpath'));
ProjDir  = fileparts(fileparts(CurrDir));
InputDir = fullfile(ProjDir, 'Input');
FlowDir = fileparts(CurrDir);
InitDir = fullfile(FlowDir, 'Steady', 'Output', 'results.mat');

%% Load Steady state Result
pres_init = load(InitDir);
steady_pres  = pres_init.output(end).pressure;

%% Import Grid
t_grid = tic;

grid = Grid(); % To initialize a mesh object
mesh_file = fullfile(InputDir, 'mesh.vtk');
grid.importMesh(mesh_file); % To import the actual mesh

%% Find Boundary nodes
bottomNodes = unique(grid.getSurfNodes(find(grid.surfaces.tag == 2)));
sideNodesAll = unique(grid.getSurfNodes(find(grid.surfaces.tag == 1)));
sideNodes = setdiff(sideNodesAll, bottomNodes);

%% Set Dirichlet Values
% We assume that the Piezometric heads 
% for the nodes on the boundary to be 0

flow_init = -grid.coordinates(:, 3)*9.81;
flow_init = flow_init([sideNodes;bottomNodes]);


%% Import Material Properties
material_file = fullfile(InputDir, 'materials.xml');
mat = Materials(material_file);
t1 = toc(t_grid);
fprintf('Grid input %.2f s\n', t1)
%% Simulation Parameters
t_simpara = tic;
simPara_file = fullfile(InputDir, 'simparam_trans.xml');
simulationInput = readInput(simPara_file);
simulationInput.LinearSolver = struct('useMatlab', 0);
simParam = SimulationParameters(simulationInput);
t2 = toc(t_simpara);
fprintf('Simulation parameters input %.2f s\n', t2);

%% Dirichlet Boundary Conditions
t_bc_flow_dir = tic;   

bc = Boundaries(grid);

bc.addBC('name'           , "Pressure_Boundary", ...
            'type'           , "Dirichlet", ...
            'field'          , "node", ...
            'variable'       , "pressure", ...
            'entityListType' , "bcList", ...
            'entityList'     , [sideNodes; bottomNodes] ...
            );

bc.addBCEvent("Pressure_Boundary", ...
                'time' , 0., ...
                'value', flow_init);

t3 = toc(t_bc_flow_dir);
fprintf('Flow DirBC imposed, %.2f s\n', t3)

%% Neumann Boundary Conditions
% The Neum. B.C. under FEM is used as setting the flux on the nodes. 

% Set up the nodes where we have extraction
INJECT_NODE = [18170, 19508];
EXTRACT_NODE = [21371, 23550, 21993];


bc.addBC( ...
    'name', "Extraction", ...
    'field', "node", ...
    'type', "Source", ...
    'variable', "pressure", ...
    'entityListType', 'bcList', ...
    'entityList', EXTRACT_NODE ...
    )
bc.addBC( ...
    'name', "Injection", ...
    'field', "node", ...
    'type', "Source", ...
    'variable', "pressure", ...
    'entityListType', 'bcList', ...
    'entityList', INJECT_NODE ...
    )

for i = eventTimes
bc.addBCEvent( ...
            "Extraction", ...
            'time', i ,...
            'value', -1e4 ...
            )

bc.addBCEvent( ...
    "Injection", ...
    'time', i, ...
    'value', 5e3)
end
% Set up the nodes where we have injection


%% Solver Settings
t_sol = tic;
domain = Discretizer('Boundaries', bc, ...
                     'Materials', mat, ...
                     'Grid', grid);

domain.addPhysicsSolver("SinglePhaseFlowFEM")
% FVTPFA -> pressure on cells;
% FEM -> pressure on nodes;

state0 = domain.getState();
state0.pressure = steady_pres;
domain.setState(state0);

out = OutState('printTimes',0:print_dt:TIME_STEPS+1,'outputFile',"Output/resultsLive",'matFileName',"Output/results");

solver = NonLinearImplicit('simulationparameters',simParam,...
    'domains',domain,...
    'output',out);

solver.simulationLoop();
t_end = toc(t_sol);
fprintf('Elapsed time in solver: %.2f s', t_end)
