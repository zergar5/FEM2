using DirectProblem.TwoDimensional.Assembling.MatrixTemplates;
using FEM2.Calculus;
using FEM2.Core.GridComponents;
using FEM2.SLAE.Preconditions;
using FEM2.SLAE.Solvers;
using FEM2.TwoDimensional.Assembling;
using FEM2.TwoDimensional.Assembling.Boundary;
using FEM2.TwoDimensional.Assembling.Global;
using FEM2.TwoDimensional.Assembling.Local;
using FEM2.TwoDimensional.Assembling.MatrixTemplates;
using FEM2.TwoDimensional.Parameters;
using System.Globalization;
using System.Security.Cryptography;
using DirectProblem.TwoDimensional;
using FEM2.Core;
using FEM2.FEM;
using FEM2.GridGenerator;
using FEM2.IO;
using FEM2.TwoDimensional;

Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

var gridI = new GridIO("../FEM2/Files/");

var nodes = gridI.ReadNodes("rz.dat");
var materialsIds = gridI.ReadMaterials("nvkat2d.dat");
var elements = gridI.ReadElements(nodes, materialsIds, "nvtr.dat");
var firstBoundaries = gridI.ReadFirstBoundaries("l1.dat");

var gridBuilder = new GridBuilder2D();

var grid = gridBuilder
    .SetNodes(nodes)
    .SetElements(elements)
    .Build();

var ironMu = 1000d * 4 * Math.PI * 1e-7;
var mu = 1d * 4 * Math.PI * 1e-7;

var materialRepository = new MaterialRepository
(
    [ironMu, mu, mu, mu],
    [0d, 0d, 1e7, -1e7]
);

var localBasisFunctionsProvider = new LocalBasisFunctionsProvider(grid, new LinearFunctionsProvider());

var localAssembler = new LocalAssembler(grid, new LocalMatrixAssembler(), materialRepository);

var globalAssembler = 
    new GlobalAssembler<Node2D>(grid, 
        new MatrixPortraitBuilder(), localAssembler, new Inserter(), new GaussExcluder());

var firstBoundaryProvider = new FirstBoundaryProvider(grid);

var firstBoundariesValues = 
    firstBoundaryProvider.GetConditions(firstBoundaries);

var solver = new MCG(new LLTPreconditioner(), new LLTSparse());

var equation = globalAssembler
    .AssembleEquation(grid)
    .ApplyFirstBoundaries(firstBoundariesValues)
    .BuildEquation();

var preconditionMatrix = globalAssembler.AllocatePreconditionMatrix();

var solution = solver
    .SetPrecondition(preconditionMatrix)
    .Solve(equation);

var femSolution = new FEMSolution(grid, solution, localBasisFunctionsProvider, new DerivativeCalculator());

Span<double> resultsAz =
[
    femSolution.CalculateAz(new Node2D(-1e-2, 1.3e-3)),
    femSolution.CalculateAz(new Node2D(-3e-3, 2.5e-3)),
    femSolution.CalculateAz(new Node2D(0, 1.6e-3)),
    femSolution.CalculateAz(new Node2D(3e-3, 2.5e-3)),
    femSolution.CalculateAz(new Node2D(9.5e-3, 2.4e-3)),
];

Span<double> resultsB =
[
    femSolution.CalculateB(new Node2D(-1e-2, 1.3e-3)),
    femSolution.CalculateB(new Node2D(-3e-3, 2.5e-3)),
    femSolution.CalculateB(new Node2D(0, 1.6e-3)),
    femSolution.CalculateB(new Node2D(3e-3, 2.5e-3)),
    femSolution.CalculateB(new Node2D(9.5e-3, 2.4e-3)),
];

foreach (var result in resultsAz)
{
    Console.WriteLine(result);
}

foreach (var result in resultsB)
{
    Console.WriteLine(result);
}