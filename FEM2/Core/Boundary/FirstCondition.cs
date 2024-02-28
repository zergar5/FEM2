using FEM2.Core.Local;

namespace FEM2.Core.Boundary;

public record struct FirstCondition(int ElementIndex, Bound Bound);
public record struct FirstConditionValue(LocalVector Values);