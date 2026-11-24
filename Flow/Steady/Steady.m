clear;
clc;

%% About the units
% m day kPa

TIME_STEPS = 1;
dt = 1;
print_dt = 1;

%% LOCATE INPUT PATH
CurrDir  = fileparts(mfilename("fullpath"));
ProjDir  = fileparts(fileparts(CurrDir));
InputDir = fullfile(ProjDir, 'Input');

%% Import Grid
t_grid = tic;

grid = Grid();               % To initialize a mesh object
mesh_file = fullfile(InputDir, 'mesh.vtk');
grid.importMesh(mesh_file); % To import the actual mesh

%% Find Boundary nodes
bottomNodes  = unique(grid.getSurfNodes(find(grid.surfaces.tag == 2)));
sideNodesAll = unique(grid.getSurfNodes(find(grid.surfaces.tag == 1)));
sideNodes    = setdiff(sideNodesAll, bottomNodes);

%% Set Dirichlet Values
flow_init = -grid.coordinates(:, 3)*9.81;
flow_init = flow_init(sideNodes);

%% Import Material Properties
material_file = fullfile(InputDir, 'materials.xml');
mat = Materials(material_file);
t1 = toc(t_grid);
fprintf('Grid input %.2f s\n', t1)
%% Simulation Parameters
t_simpara = tic;
simPara_file = fullfile(InputDir, 'simparam_steady.xml');
simulationInput = readInput(simPara_file);
simulationInput.LinearSolver = struct('useMatlab', 0);
simParam = SimulationParameters(simulationInput);
t2 = toc(t_simpara);
fprintf('Simulation parameters input %.2f s\n', t2);
%% Boundary conditions
t_bc_flow_dir = tic;   

bound = Boundaries(grid);

bound.addBC('name'           , "Pressure_Boundary", ...
            'type'           , "Dirichlet", ...
            'field'          , "node", ...
            'variable'       , "pressure", ...
            'entityListType' , "bcList", ...
            'entityList'     , sideNodes ...
            );

bound.addBCEvent("Pressure_Boundary", ...
                'time' , 0., ...
                'value', flow_init);

t3 = toc(t_bc_flow_dir);
fprintf('Flow DirBC imposed, %.2f s\n', t3)


%% Solver Settings
t_sol = tic;
domain = Discretizer('Boundaries', bound, ...
                     'Materials', mat, ...
                     'Grid', grid);

domain.addPhysicsSolver("SinglePhaseFlowFEM", "steadyState", 1)

out = OutState( ...
    'printTimes',0:print_dt:TIME_STEPS+1, ...
    'outputFile',"Output/resultsLive", ...
    'matFileName',"Output/results");

solver = NonLinearImplicit( ...
    'simulationparameters',simParam,...
    'domains',domain,...
    'output',out);


solver.simulationLoop();
t_end = toc(t_sol);
fprintf('Elapsed time in solver: %.2f s', t_end)
