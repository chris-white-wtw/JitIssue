using System.Diagnostics.CodeAnalysis;

namespace JitIssue
{
  /// <summary>
  ///   Gets the mathematical constants used across files.
  /// </summary>
  /// <value>
  ///   The constants used across files.
  /// </value>
  /// <remarks>
  ///   Constants used in only one file are better placed in that file.<seealso cref="DateTimeConstants" />
  /// </remarks>
  [SuppressMessage("ReSharper", "InconsistentNaming")]
  internal static class Constants
  {
    /// <summary>
    ///   Gets the machine precision.
    /// </summary>
    /// <remarks>
    ///   2**-53 from Delphi.
    /// </remarks>
    internal const double MachinePrecision = 1.11022302462515654042e-16;

    /// <summary>
    ///   Gets the maximum log result.
    /// </summary>
    /// <value>
    ///   MAXLOG and log(MAXNUM) from Delphi.
    /// </value>
    internal const double MaximumLogResults = 7.09782712893383996732e2;

    /// <summary>
    ///   Gets the minimum log result
    /// </summary>
    /// <value>
    ///   MINLOG from Delphi (log(2**-1075)).
    /// </value>
    internal const double MinimumLogResults = -7.451332191019412076235e2;
  }
}