clear;
clc;

%% Paths and units: m, day, kPa, kN/m^3
CurrDir = fileparts(mfilename('fullpath'));
transientDir = fullfile(fileparts(CurrDir), 'Flow', 'Transient');

ProjDir  = fileparts(CurrDir);
InputDir = fullfile(ProjDir, 'Input');

%% Input Material Properties and Simulation Parameters Files
materialFile = fullfile(InputDir, 'materials.xml');
simulationInput = readInput(fullfile(InputDir, 'simparam_trans.xml'));
simulationInput.LinearSolver = struct('useMatlab', 0);
simParam = SimulationParameters(simulationInput);
outputDir = fullfile(CurrDir, 'Output');

%% Read the flow mesh and its nodal pressure history
grid = Grid();
grid.importMesh(fullfile(InputDir, 'mesh.vtk'));
grid.processGeometry();
flow = load(fullfile(transientDir, 'Output', 'results.mat'), 'output');
[pressureTimes, timeOrder] = sort([flow.output.time]);
pressureNodes = [flow.output(timeOrder).pressure];

% For this four-node tetrahedral mesh, the nodal mean is the exact
% element-average pressure of the linear FEM pressure field.
cellNodes = grid.getCellNodes((1:grid.cells.num)');
pressureCells = zeros(grid.cells.num, numel(pressureTimes));

for k = 1:numel(pressureTimes)
    nodalPressure = pressureNodes(:, k);
    pressureCells(:, k) = mean(nodalPressure(cellNodes), 2);
end

p0Cell = pressureCells(:, 1);
deltaPressure = pressureCells - p0Cell;

%% Materials and the tagged bottom boundary
materials = Materials(materialFile);
bottomNodes = unique(grid.getSurfNodes(find(grid.surfaces.tag == 2)));
sideNodesAll = unique(grid.getSurfNodes(find(grid.surfaces.tag == 1)));
sideNodes = setdiff(sideNodesAll, bottomNodes);

%% Initial effective stress under gravity and the full initial pore pressure
fprintf('Computing initial effective stress...\n');
[sigma0, equilibriumState, gpMap] = computeInitialEffectiveStress( ...
    grid, materialFile, p0Cell, bottomNodes, simulationInput);

%% Mechanical boundary conditions and one-way pressure loading
bc = Boundaries(grid);
bc.addBC( ...
    'name', "FixedBottom", ...
    'type', "Dirichlet", ...
    'field', "node", ...
    'variable', "displacements", ...
    'entityListType', "bcList", ...
    'entityList', bottomNodes, ...
    'components', ["x", "y", "z"]);

bc.addBCEvent("FixedBottom", 'time', simParam.tIni, 'value', 0.);

% Side and top boundaries have no prescribed displacement or added traction.
% Poromechanics.volumeForce integrates the supplied scalar against B^T I;
% therefore the supplied pressure load includes the material Biot coefficient.
bc.addBC( ...
    'name', "PressureChange", ...
    'type', "volumeForce", ...
    'field', "cell", ...
    'variable', "displacements", ...
    'entityListType', "bcList", ...
    'entityList', (1:grid.cells.num)');

for k = 1:numel(pressureTimes)
    bc.addBCEvent("PressureChange", 'time', pressureTimes(k), ...
        'value', deltaPressure(:, k));
end

%% Displacements are measured from the initial equilibrium configuration
domain = Discretizer('Boundaries', bc, 'Materials', materials, 'Grid', grid);
domain.addPhysicsSolver("Poromechanics");
state0 = domain.getState();
state0.stress = sigma0;
state0.displacements(:) = 0.;
state0.strain(:) = 0.;
domain.setState(state0);

printTimes = unique([simParam.tIni, ...
    pressureTimes(pressureTimes >= simParam.tIni & pressureTimes <= simParam.tMax), ...
    simParam.tMax]);

out = OutState('printTimes', printTimes, 'solvePrintTimes', 1, ...
    'outputFile', fullfile(outputDir, 'resultsLive'), ...
    'matFileName', fullfile(outputDir, 'results'));

solver = NonLinearImplicit('simulationparameters', simParam, ...
    'domains', domain, 'output', out);
solver.simulationLoop();
save(fullfile(outputDir, 'initialEffectiveStress.mat'), ...
    'sigma0', 'equilibriumState', 'gpMap', 'p0Cell', 'bottomNodes');
fprintf('Geomechanical results saved in %s\n', outputDir);
