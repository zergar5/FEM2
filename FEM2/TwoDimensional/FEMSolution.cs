using FEM2.Calculus;
using FEM2.Core;
using FEM2.Core.GridComponents;
using FEM2.FEM;
using FEM2.SLAE;
using FEM2.TwoDimensional.Assembling.Local;
using Vector = FEM2.Core.Base.Vector;

namespace FEM2.TwoDimensional;

public class FEMSolution
{
    private readonly Grid<Node2D> _grid;
    private readonly Vector _solution;
    private readonly LocalBasisFunctionsProvider _localBasisFunctionsProvider;
    private readonly DerivativeCalculator _derivativeCalculator;

    public FEMSolution
    (
        Grid<Node2D> grid,
        Vector solution,
        LocalBasisFunctionsProvider localBasisFunctionsProvider,
        DerivativeCalculator derivativeCalculator
    )
    {
        _grid = grid;
        _solution = solution;
        _localBasisFunctionsProvider = localBasisFunctionsProvider;
        _derivativeCalculator = derivativeCalculator;
    }

    public double CalculateAz(Node2D point)
    {
        if (AreaHas(point))
        {
            var element = _grid.Elements.First(x => ElementHas(x, point));

            var basisFunctions = _localBasisFunctionsProvider.GetBilinearFunctions(element);

            var sum = 0d;

            sum += element.NodesIndexes
                .Select((t, i) => _solution[t] * basisFunctions[i].Calculate(point))
                .Sum();

            CourseHolder.WriteAz(point, sum);

            return sum;
        }

        CourseHolder.WriteAreaInfo();
        CourseHolder.WriteAz(point, 0);
        return 0;
    }

    public double CalculateB(Node2D point)
    {
        if (AreaHas(point))
        {
            var element = _grid.Elements.First(x => ElementHas(x, point));

            var basisFunctions = _localBasisFunctionsProvider.GetBilinearFunctions(element);

            var sumX = 0d;
            var sumY = 0d;

            sumX += element.NodesIndexes
                .Select((t, i) => 
                    _solution[t] * _derivativeCalculator.Calculate(basisFunctions[i], point, 'x'))
                .Sum();

            sumY += element.NodesIndexes
                .Select((t, i) =>
                    _solution[t] * _derivativeCalculator.Calculate(basisFunctions[i], point, 'y'))
                .Sum();

            var b = Math.Sqrt(sumX * sumX + sumY * sumY);

            CourseHolder.WriteAz(point, b);

            return b;
        }

        CourseHolder.WriteAreaInfo();
        CourseHolder.WriteAz(point, 0);
        return 0;
    }

    private bool ElementHas(Element element, Node2D node)
    {
        var lowerLeftCorner = _grid.Nodes[element.NodesIndexes[0]];
        var upperRightCorner = _grid.Nodes[element.NodesIndexes[^1]];
        return (node.X > lowerLeftCorner.X ||
                Math.Abs(node.X - lowerLeftCorner.X) < MethodsConfig.EpsDouble) &&
               (node.Y > lowerLeftCorner.Y ||
                Math.Abs(node.Y - lowerLeftCorner.Y) < MethodsConfig.EpsDouble) &&
               (node.X < upperRightCorner.X ||
                Math.Abs(node.X - upperRightCorner.X) < MethodsConfig.EpsDouble) &&
               (node.Y < upperRightCorner.Y ||
                Math.Abs(node.Y - upperRightCorner.Y) < MethodsConfig.EpsDouble);
    }

    private bool AreaHas(Node2D node)
    {
        var lowerLeftCorner = _grid.Nodes[0];
        var upperRightCorner = _grid.Nodes[^1];
        return (node.X > lowerLeftCorner.X ||
                Math.Abs(node.X - lowerLeftCorner.X) < MethodsConfig.EpsDouble) &&
               (node.Y > lowerLeftCorner.Y ||
                Math.Abs(node.Y - lowerLeftCorner.Y) < MethodsConfig.EpsDouble) &&
               (node.X < upperRightCorner.X ||
                Math.Abs(node.X - upperRightCorner.X) < MethodsConfig.EpsDouble) &&
               (node.Y < upperRightCorner.Y ||
                Math.Abs(node.Y - upperRightCorner.Y) < MethodsConfig.EpsDouble);
    }
}