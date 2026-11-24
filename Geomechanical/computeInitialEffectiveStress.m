function [sigma0, equilibriumState, gpMap] = computeInitialEffectiveStress( ...
        grid, materialFile, p0Cell, bottomNodes, simulationInput)
% Compute the initial effective stress in a separate static equilibrium solve.
% specificWeight in PorousRock is interpreted as total bulk unit weight.
% Initial loading displacements are returned for reference only; the main
% simulation measures subsequent displacement from this equilibrium state.

input = readInput(materialFile);
for k = 1:numel(input.Solid)
    lawNames = fieldnames(input.Solid(k).Constitutive);
    law = input.Solid(k).Constitutive.(lawNames{1});
    input.Solid(k).Constitutive = struct('Elastic', struct( ...
        'youngModulus', law.youngModulus, 'poissonRatio', law.poissonRatio));
    % PreStressMechanics reads Solid.specificWeight for gravity.
    input.Solid(k).specificWeight = input.Solid(k).PorousRock.specificWeight;
end
initialMaterials = Materials(input);

bc = Boundaries(grid);
bc.addBC('name', "InitialFixedBottom", 'type', "Dirichlet", ...
    'field', "node", 'variable', "displacements", ...
    'entityListType', "bcList", 'entityList', bottomNodes, ...
    'components', ["x", "y", "z"]);
bc.addBCEvent("InitialFixedBottom", 'time', 0., 'value', 0.);

bc.addBC('name', "InitialPorePressure", 'type', "volumeForce", ...
    'field', "cell", 'variable', "displacements", ...
    'entityListType', "bcList", 'entityList', (1:grid.cells.num)');
bc.addBCEvent("InitialPorePressure", 'time', 0., ...
    'value', p0Cell);

initialDomain = Discretizer('Boundaries', bc, ...
    'Materials', initialMaterials, 'Grid', grid);
% Total gravity and the full initial pore pressure are applied separately;
% do not subtract hydrostatic buoyancy a second time.
initialDomain.addPhysicsSolver("PreStressMechanics", 'hydrostaticPressure', 0);

% Auxiliary loading time is independent of the physical simulation duration.
initialInput = simulationInput;
initialInput.Start = 0.;
initialInput.End = 1.;
initialInput.DtInit = 1.;
initialInput.DtMin = 1.;
initialInput.DtMax = 1.;
initialInput.incrementFactor = 1.;
initialInput.LinearSolver = struct('useMatlab', 0);
initialSolver = NonLinearImplicit( ...
    'simulationparameters', SimulationParameters(initialInput), ...
    'domains', initialDomain, 'output', OutState());
initialSolver.simulationLoop();

equilibriumState = initialDomain.getState();
sigma0 = equilibriumState.stress;
gpMap = initialDomain.gpMap;
end
