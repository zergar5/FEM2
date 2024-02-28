using FEM2.Core.GridComponents;

namespace FEM2.FEM;

public class CourseHolder
{
    public static void GetInfo(int iteration, double residual)
    {
        Console.Write($"Iteration: {iteration}, residual: {residual:E14}                                   \r");
    }

    public static void WriteAz(Node2D point, double value)
    {
        Console.WriteLine($"({point.X:F4}, {point.Y:F4}) Az = {value}");
    }

    public static void WriteB(Node2D point, double value)
    {
        Console.WriteLine($"({point.X:F4}, {point.Y:F4}) B = {value}");
    }

    public static void WriteAreaInfo()
    {
        Console.WriteLine("Point not in area or time interval");
    }
}