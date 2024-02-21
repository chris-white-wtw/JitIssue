namespace JitIssue;

internal class Program
{
  record struct IibiParams(double aa, double bb, double tolerance, double yy0);

  static void Main(string[] args)
  {
    var random = new Random(1);

    int i = 0;
    IibiParams p = default;
    try
    {
      for (i = 0; i < 300000; i++)
      {
        p.aa = 1;
        p.bb = 2;
        p.tolerance = 0;
        p.yy0 = random.NextDouble();

        CephesPort.InverseIncompleteBetaIntegral(p.aa, p.bb, new(p.tolerance), p.yy0);
      }
    }
    catch (Exception ex)
    {
      Console.Error.WriteLine($"Exception on iteration {i}, params {p}");
      Console.Error.WriteLine($"aplusb={CephesPort.lastAPlusB}, a={CephesPort.lastA}, b={CephesPort.lastB}, a for power series={CephesPort.lastAForPowerSeries}");
      Console.Error.WriteLine(ex);
    }
  }
}