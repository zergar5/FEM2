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
using DirectProblem.FEM;
using DirectProblem.TwoDimensional;
using FEM2.Core;
using FEM2.GridGenerator;
using FEM2.IO;
using FEM2.TwoDimensional;

Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

var gridI = new GridIO("../FEM2/Files/");

var nodes = gridI.ReadNodes("rz.dat");
var materialsIds = gridI.ReadMaterials("nvkat2d.dat");
var elements = gridI.ReadElements(nodes, materialsIds, "nvtr.dat");

var gridBuilder = new GridBuilder2D();

var grid = gridBuilder
    .SetNodes(nodes)
    .SetElements(elements)
    .Build();

var materialRepository = new MaterialRepository
(
    [1000d * 4 * Math.PI * 1e-7, 1d * 4 * Math.PI * 1e-7, 1d * 4 * Math.PI * 1e-7, 1d * 4 * Math.PI * 1e-7],
    [0d, 0d, 1e7, -1e7]
);

var localMatrixAssembler = new LocalMatrixAssembler();

var localBasisFunctionsProvider = new LocalBasisFunctionsProvider(grid, new LinearFunctionsProvider());

var derivativeCalculator = new DerivativeCalculator();

var localAssembler = new LocalAssembler(grid, localMatrixAssembler, materialRepository);

var inserter = new Inserter();
var globalAssembler = new GlobalAssembler<Node2D>(grid, new MatrixPortraitBuilder(), localAssembler, inserter, new GaussExcluder());

var firstBoundaryProvider = new FirstBoundaryProvider(grid);

var firstBoundaries = firstBoundaryProvider.GetConditions(firstBoundaryProvider.GetArrays(11, 6));

var lltPreconditioner = new LLTPreconditioner();
var solver = new MCG(lltPreconditioner, new LLTSparse());

var equation = globalAssembler
    .AssembleEquation(grid)
    .ApplyFirstConditions(firstBoundaries)
    .BuildEquation();

var preconditionMatrix = globalAssembler.AllocatePreconditionMatrix();

var solution = solver
    .SetPrecondition(preconditionMatrix)
    .Solve(equation);

var femSolution = new FEMSolution(grid, solution, localBasisFunctionsProvider);

Span<double> results =
[
    femSolution.CalculateAz(new Node2D(-1e-2, 1.3e-3)),
    femSolution.CalculateAz(new Node2D(-3e-3, 2.5e-3)),
    femSolution.CalculateAz(new Node2D(0, 1.6e-3)),
    femSolution.CalculateAz(new Node2D(3e-3, 2.5e-3)),
    femSolution.CalculateAz(new Node2D(9.5e-3, 2.4e-3)),
];

foreach (var result in results)
{
    Console.WriteLine(result);
}