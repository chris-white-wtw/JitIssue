using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.Runtime.CompilerServices;

namespace JitIssue
{
  /// <summary>
  ///   Methods adapted from those in the Cephes Math Library.
  /// </summary>
  /// <remarks>
  ///   <para>Adapted from Cephes Math Library Release 2.9:  November, 2000</para>
  ///   <para>Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier</para>
  ///   <para>And since placed in the public domain</para>
  /// </remarks>
  [SuppressMessage("Microsoft.StyleCop.CSharp.DocumentationRules", "SA1600:ElementsMustBeDocumented", Justification = "There are lots of fields here that don't need documenting")]
  [SuppressMessage("ReSharper", "InconsistentNaming", Justification = "Ported code")]
  internal static class CephesPort
  {
    private const double cBig = 4.503599627370496e15;

    private const double cBigInv = 2.22044604925031308085e-16; // reciprocal of cBig

    private const double cRootMachinePrecision = 1.490116119384765625E-8;

    private const double cMaxGam = 171.624376956302725;

    public struct Tolerance
    {
      public readonly double Value;

      public Tolerance(double val)
      {
        Value = val;
      }

    }
    /// <summary>
    /// The value to pass to methods taking a tolerance if you want the method to use a sensible default
    /// </summary>
    public static readonly Tolerance DefaultTolerance = new Tolerance(0);

    /// <summary>
    ///   Stirling's Formula expansion of log gamma.
    /// </summary>
    private static readonly double[] cA =
    {
      8.11614167470508450300E-4,
      -5.95061904284301438324E-4,
      7.93650340457716943945E-4,
      -2.77777777730099687205E-3,
      8.33333333333331927722E-2
    };

    /// <summary>
    ///   Log gamma at 2.
    /// </summary>
    private static readonly double[] cB =
    {
      -1.37825152569120859100E3,
      -3.88016315134637840924E4,
      -3.31612992738871184744E5,
      -1.16237097492762307383E6,
      -1.72173700820839662146E6,
      -8.53555664245765465627E5
    };

    /// <summary>
    ///   Log gamma at 3.
    /// </summary>
    private static readonly double[] cC =
    {
      -3.51815701436523470549E2,
      -1.70642106651881159223E4,
      -2.20528590553854454839E5,
      -1.13933444367982507207E6,
      -2.53252307177582951285E6,
      -2.01889141433532773231E6
    };

    private static readonly double[] cPgam =
    {
      1.60119522476751861407E-4,
      1.19135147006586384913E-3,
      1.04213797561761569935E-2,
      4.76367800457137231464E-2,
      2.07448227648435975150E-1,
      4.94214826801497100753E-1,
      9.99999999999999996796E-1
    };

    private static readonly double[] cQgam =
    {
      -2.31581873324120129819E-5,
      5.39605580493303397842E-4,
      -4.45641913851797240494E-3,
      1.18139785222060435552E-2,
      3.58236398605498653373E-2,
      -2.34591795718243348568E-1,
      7.14304917030273074085E-2,
      1.00000000000000000320E0
    };

    private static readonly double[] cStirling =
    {
      7.87311395793093628397E-4,
      -2.29549961613378126380E-4,
      -2.68132617805781232825E-3,
      3.47222221605458667310E-3,
      8.33333333333482257126E-2
    };

    /// <summary>
    ///   Calculates the left tail of the incomplete gamma function:
    /// </summary>
    /// <param name="a">The a parameter.</param>
    /// <param name="x">The x parameter.</param>
    /// <returns>
    ///   The left tail of the incomplete gamma function with the supplied parameters.
    /// </returns>
    /// <remarks>This is the Cephes <c>igam</c> method.</remarks>
    internal static double IncompleteGamma(double a, double x)
    {
      if (x <= 0 || a <= 0)
      {
        return 0;
      }

      if (1 < x && a < x)
      {
        return 1 - ComplementedIncompleteGamma(a, x);
      }

      // Compute x**a * exp(-x) / gamma(a)
      double ax = (a * Math.Log(x)) - x - LogGamma(a);
      if (ax < -Constants.MaximumLogResults)
      {
        return 0;
      }

      ax = Math.Exp(ax);

      // power series
      double r = a;
      double c = 1;
      double ans = 1;
      double oldR;

      double machinePrecision = Constants.MachinePrecision; // we're going to call this a lot, so cache it

      do
      {
        oldR = r;
        r += 1;

        if (r == oldR)
        {
          // because r was large (> 9e15), adding 1 made no difference
          // if we don't break, the outer loop could go round and round for a very long time
          // and we won't be able to trust the result anyway
          throw new ArithmeticException("Overflow calculating power series");
        }

        c *= x / r;
        ans += c;
      } while (ans * machinePrecision < c);  // ans and c both +ve, so this is equivalent to and faster than mp<c/ans

      double result = ans * ax / a;

      if (result < 0 || 1 < result)
      {
        throw new ArithmeticException("Overflow calculating power series (2)");
      }

      return result;
    }

    /// <summary>
    ///   Calculates the inverse of the complemented incomplete Gamma integral.
    /// </summary>
    /// <param name="a">The a parameter.</param>
    /// <param name="y0">The y0 parameter.</param>
    /// <returns>
    ///   <para>The inverse of the complemented incomplete Gamma function.</para>
    ///   <para>
    ///     Returns <c>x</c> such that <see cref="ComplementedIncompleteGamma" />(<paramref name="a" />, x) =
    ///     <paramref name="y0" />.
    ///   </para>
    /// </returns>
    /// <remarks>
    ///   This is the Cephes <c>igami</c> method.
    /// </remarks>
    internal static double InverseComplementedIncompleteGamma(double a, double y0)
    {
      // TODO (Fix) I believe this function is only valid in the right-hand tail (i.e. where y0 <= 0.5)
      //    but it's being called with 0 <= y0 <= 1. This needs checking, with reference to
      //    http://search.cpan.org/~rkobes/Math-Cephes-0.45/lib/Math/Cephes.pod and
      //    http://www.netlib.org/cephes/doubldoc.html#igami
      //    VSTS 83949

      // bound the solution
      double x0 = double.MaxValue;
      double yl = 0; // low bound of y
      double x1 = 0;
      double yh = 1; // high bound of y
      double dithresh = 5 * Constants.MachinePrecision;

      // approximation to inverse function
      double d = 1 / (9 * a);
      double y = 1 - d - (MorosInvNormal(y0) * Math.Sqrt(d)); // was ndtri
      double x = a * y * y * y;

      double lgm = LogGamma(a);

      for (int j = 0; j <= 9; j++)
      {
        if (x0 < x || x < x1)
        {
          goto ihalve;
        }

        y = ComplementedIncompleteGamma(a, x);

        if ((y < yl) || (yh < y))
        {
          goto ihalve;
        }

        if (y < y0)
        {
          x0 = x;
          yl = y;
        }
        else
        {
          x1 = x;
          yh = y;
        }

        // compute the derivative of the function at this point
        d = ((a - 1.0) * Math.Log(x)) - x - lgm;
        if (d < -Constants.MaximumLogResults)
        {
          goto ihalve;
        }

        d = -Math.Exp(d);

        // compute the step to the next approximation of x
        d = (y - y0) / d;
        if (CheckRatio(d, x, Constants.MachinePrecision))
        {
          return x;
        }

        x -= d;
      }

    //! Resort to interval halving if Newton iteration did not converge.
    ihalve:

      d = 0.0625;
      if (x0 == double.MaxValue)
      {
        if (x <= 0)
        {
          x = 1;
        }

        while (true /* was x0 == double.MaxValue , but x0 is only changed just before a break statement*/)
        {
          x = (1 + d) * x;
          y = ComplementedIncompleteGamma(a, x);
          if (y < y0)
          {
            x0 = x;
            yl = y;
            break;
          }

          d = d + d;

          if (double.IsPositiveInfinity(d))
          {
            throw new ArithmeticException("Arithmetic overflow");
          }
        }
      }

      d = 0.5;
      int dir = 0;

      int i = 0;
      for (; i < 400; i++)
      {
        x = x1 + (d * (x0 - x1));
        y = ComplementedIncompleteGamma(a, x);
        if (x <= 0)
        {
          break;
        }

        //lgm = (x0 - x1) / (x1 + x0);
        //if (Math.Abs(lgm) < dithresh)
        if (CheckRatio(x0 - x1, x1 + x0, dithresh))
        {
          break;
        }

        //lgm = (y - y0) / y0;
        if (CheckRatio(y - y0, y0, dithresh))
        {
          break;
        }

        if (y0 <= y)
        {
          x1 = x;
          yh = y;
          if (dir < 0)
          {
            dir = 0;
            d = 0.5;
          }
          else if (1 < dir)
          {
            d = (0.5 * d) + 0.5;
          }
          else
          {
            d = (y0 - yl) / (yh - yl);
          }

          dir = dir + 1;
        }
        else
        {
          x0 = x;
          yl = y;
          if (0 < dir)
          {
            dir = 0;
            d = 0.5;
          }
          else if (dir < -1)
          {
            d = 0.5 * d;
          }
          else
          {
            d = (y0 - yl) / (yh - yl);
          }

          dir = dir - 1;
        }
      }

      return x;
    }



    private static class InverseIncompleteBetaIntegralLargeApproximation
    {
      private static readonly double[] cEps1C1 = { 75, 202, 188, -888, -1345, 118, 138 };
      private static readonly double[] cEps1C2 = { 7, 21, 70, 26, -93, -31 };
      private static readonly double[] cEps1C3 = { 1, -13, 69, 167, 46 };
      private static readonly double[] cEps1C4 = { 1, 9, 21, 5 };
      private static readonly double[] cEps2C1 = { 11053, 53308, 117010, 163924, 116188, -258428, -677042, -481940, -105497 };
      private static readonly double[] cEps2C2 = { 2132, 7915, 16821, 35066, 87490, 141183, 95993, 21640 };
      private static readonly double[] cEps2C3 = { 35, -154, -623, -1636, -3983, -3514, -925 };
      private static readonly double[] cEps2C4 = { 28, 131, 402, 581, 208 };
      private static readonly double[] cEps3C1 = { 116932, 819281, 2378172, 4341330, 6806004, 10622748, 18739500, 30651894, 30869976, 15431867, 2919016 };
      private static readonly double[] cEps3C2 = { 442043, 2054169, 3803094, 3470754, 2141568, -2393568, -19904934, -34714674, 23128299, -5253353 };
      private static readonly double[] cEps3C3 = { 3592, 8375, -1323, -29198, -89578, -154413, -116063, -29632 };

      private static double Refine(double x, double p, double a, double b, Tolerance tolerance)
      {
        double p0 = IncompleteBetaIntegral(a, b, tolerance, x);

        double result = x;

        if ((p > 0) && (p0 > 0) && (p0 < 1) && (p != p0))
        {
          double z = -Math.Log(1 / p - 1);
          double z0 = -Math.Log(1 / p0 - 1);

          if (z != z0)
          {
            double invDeriv = Math.Exp(-(LogGamma(a + b) - LogGamma(a) - LogGamma(b) + (b - 1) * Math.Log(1 - x) + (a - 1) * Math.Log(x))) * (p0 - p0 * p0);

            result = x + invDeriv * (z - z0);
          }
        }

        result = Math.Max(Math.Min(result, 1.0), 0.0);
        return result;
      }
      internal static double DoApproximation(double aa, double bb, Tolerance tolerance, double yy0)
      {
        double a;
        double b;
        // From 'Switch?' base unit:
        // Approximation is valid only for a/b>>1. If b>a then make use of the identiy Beta(a,b,p)=1-Beta(b,a,1-p)
        bool reverse = (bb > aa);
        double p;
        if (reverse)
        {
          a = bb;
          b = aa;
          p = 1 - yy0;
        }
        else
        {
          a = aa;
          b = bb;
          p = yy0;
        }


        // From 'InverseBeta Approximation' base unit:
        // Uses asymptotic approximation for inverse of incomplete beta function for
        // large a and a>>b from Temme, Asymptotic Inversion of the Incomplete Beta Function,
        // Journal of Computational and Applied Mathematics, 41 (1992) 145-157.
        // Terms of order (b/a) are kept.
        // If a<b then the approximation is for large b and b>>a
        // NB: Solver used in the model was not required
        double eta0 = Math.Max(1 / a * InverseComplementedIncompleteGamma(b, p), 1e-13);

        double mu = b / a;

        double w = Math.Sqrt(1 + mu);

        double epsilon1 = PolyEval((eta0 - mu) / (w * (w + 1)),
                            new[]
                            {
                            -PolyEval(w, cEps1C1) / 272160,
                            -PolyEval(w, cEps1C2) / 6480,
                            -PolyEval(w, cEps1C3) / 1620,
                            PolyEval(w, cEps1C4) / 36,
                            (w + 2) * (w - 1) / 3
                            }) / w;

        double epsilon2 = PolyEval((eta0 - mu) / (w * (w + 1)),
                            new[]
                            {
                            -PolyEval(w, cEps2C1) / 14696640,
                            -PolyEval(w, cEps2C2) / 816480,
                            -PolyEval(w, cEps2C3) / 12960,
                            PolyEval(w, cEps2C4) * (w - 1) / 1620
                            }) / ((w + 1) * w * w * w);

        double epsilon3 = -PolyEval((eta0 - mu) / (w * (w + 1)),
                            new[]
                            {
                            PolyEval(w, cEps3C1) / 146966400,
                            PolyEval(w, cEps3C2) / 146966400,
                            PolyEval(w, cEps3C3) * (w - 1) / 816480
                            }) / ((w + 1) * (w + 1) * w * w * w * w * w);


        double x = Math.Exp(-(eta0 + epsilon1 / (a) + epsilon2 / (a * a) + epsilon3 / (a * a * a)));

        x = Refine(x, p, a, b, tolerance);

        var result = Refine(x, p, a, b, tolerance); // a second Newton-Raphson step was necessary for accurary

        // From 'Switch Back' base unit:
        if (reverse)
        {
          result = 1.0 - result;
        }

        return result;
      }
    }

    /// <summary>
    ///   Calculates the inverse of the incomplete beta integral.
    /// </summary>
    /// <param name="aa">The alpha.</param>
    /// <param name="bb">The beta.</param>
    /// <param name="tolerance">The tolerance to pass to the iterative calculation</param>
    /// <param name="yy0">The value of y.</param>
    /// <returns>The inverse of the incomplete beta integral.</returns>
    /// <remarks>
    ///   <p>This is the Cephes <c>incbi</c> function.</p>
    ///   <p>Given y, the function finds x such that incbet(a, b, x) = y.</p>
    ///   <p>The routine performs interval halving or Newton iterations to find the root of incbet(a, b, x) - y = 0.</p>
    /// </remarks>
    internal static double InverseIncompleteBetaIntegral(double aa, double bb, Tolerance tolerance, double yy0)
    {
      if (yy0 <= 0)
      {
        return 0;
      }

      if (yy0 >= 1.0)
      {
        return 1;
      }

      double a;
      double b;
      double x;

      // Conditions as per e-mail from JamesN attached to VSTS: 23550 (InvBeta.msg)
      if ((aa > 1000) && (aa > bb * 1000) && (yy0 < 1e-5) ||
          ((bb > 1000) && (bb > aa * 1000) && ((1 - yy0) < 1e-5)))
      {
        return InverseIncompleteBetaIntegralLargeApproximation.DoApproximation(aa, bb, tolerance, yy0);
      }

      // variables used in Newton iterations and interval halving
      double x0 = 0;
      double yl = 0;
      double x1 = 1;
      double yh = 1;

      // flags
      bool hasDoneNewtonIteration = false;
      bool isReversed;

      // other variables
      double dithresh;

      double y0;
      double y;

      if (aa <= 1.0 || bb <= 1.0)
      {
        dithresh = 1.0e-6;
        isReversed = false;
        a = aa;
        b = bb;
        y0 = yy0;
        x = a / (a + b);
        y = IncompleteBetaIntegral(a, b, tolerance, x);
        goto ihalve;
      }

      dithresh = 1.0e-4;
      // approximation to inverse function

      double yp = -MorosInvNormal(yy0); //-ndtri(yy0);

      if (yy0 > 0.5)
      {
        isReversed = true;
        a = bb;
        b = aa;
        y0 = 1 - yy0;
        yp = -yp;
      }
      else
      {
        isReversed = false;
        a = aa;
        b = bb;
        y0 = yy0;
      }

      double lgm = ((yp * yp) - 3.0) / 6.0;

      x = 2.0 / ((1.0 / ((2.0 * a) - 1.0)) + (1.0 / ((2.0 * b) - 1.0)));

      double d = (yp * (Math.Sqrt(x + lgm) / x))
                 - (((1.0 / ((2.0 * b) - 1.0)) - (1.0 / ((2.0 * a) - 1.0)))
                    * (lgm + (5.0 / 6.0) - (2.0 / (3.0 * x))));

      d = 2.0 * d;

      if (d < Constants.MinimumLogResults)
      {
        x = 0;
        goto done;
        //x goto under;
      }

      x = a / (a + (b * Math.Exp(d)));

      y = IncompleteBetaIntegral(a, b, tolerance, x);

      //yp = (y - y0) / y0;

      if (CheckRatio(y - y0, y0, 0.2))
      {
        goto newt;
      }

    // Resort to interval halving if not close enough
    //!++ ihalve:
    ihalve:

      int dir = 0;
      double di = 0.5;
      for (int i = 0; i <= 99; i++)
      {

        if (i != 0)
        {
          x = x0 + (di * (x1 - x0));

          if (x == 1.0)
          {
            x = 1.0 - Constants.MachinePrecision;
            Debug.Assert(x < 1, "Should have subtracted something just large enough to make a difference");
          }

          if (x == 0.0)
          {
            di = 0.5;
            x = x0 + (di * (x1 - x0));
            if (x == 0.0)
            {
              goto done;
              //x goto under;
            }
          }

          y = IncompleteBetaIntegral(a, b, tolerance, x);
          if (CheckRatio(x1 - x0, x1 + x0, dithresh))
          {
            goto newt;
          }


          if (CheckRatio(y - y0, y0, dithresh))
          {
            goto newt;
          }
        }

        if (y < y0)
        {
          x0 = x;
          yl = y;

          if (dir < 0)
          {
            dir = 0;
            di = 0.5;
          }
          else if (dir > 3)
          {
            di = 1.0 - ((1.0 - di) * (1.0 - di));
          }
          else if (dir > 1)
          {
            di = (0.5 * di) + 0.5;
          }
          else
          {
            di = (y0 - y) / (yh - yl);
          }

          dir++;

          if (x0 > 0.75)
          {
            if (isReversed)
            {
              isReversed = false;
              a = aa;
              b = bb;
              y0 = yy0;
            }
            else
            {
              isReversed = true;
              a = bb;
              b = aa;
              y0 = 1.0 - yy0;
            }

            x = 1.0 - x;
            y = IncompleteBetaIntegral(a, b, tolerance, x);
            x0 = 0.0;
            yl = 0.0;
            x1 = 1.0;
            yh = 1.0;
            goto ihalve;
          }
        }
        else
        {
          x1 = x;

          if (isReversed && (x1 < Constants.MachinePrecision))
          {
            x = 0.0;
            goto done;
          }

          yh = y;

          if (dir > 0)
          {
            dir = 0;
            di = 0.5;
          }
          else if (dir < -3)
          {
            di = di * di;
          }
          else if (dir < -1)
          {
            di = 0.5 * di;
          }
          else
          {
            di = (y - y0) / (yh - yl);
          }

          dir--;
        }
      }

      if (x0 >= 1.0)
      {
        x = 1.0 - Constants.MachinePrecision;
        goto done;
      }

      if (x <= 0.0)
      {
        x = 0.0;
        goto done;
      }

    //!+ newt: Newton iteration
    newt:

      if (hasDoneNewtonIteration)
      {
        goto done;
      }

      double newtonTolerance = tolerance.Value == 0 ? 128 * Constants.MachinePrecision : tolerance.Value;
      hasDoneNewtonIteration = true;
      lgm = LogGamma(a + b) - LogGamma(a) - LogGamma(b);

      for (int i = 0; i <= 7; i++)
      {
        // Compute the function at this point
        if (i != 0)
        {
          y = IncompleteBetaIntegral(a, b, tolerance, x);
        }

        if (y < yl)
        {
          x = x0;
          y = yl;
        }
        else if (y > yh)
        {
          x = x1;
          y = yh;
        }
        else if (y < y0)
        {
          x0 = x;
          yl = y;
        }
        else
        {
          x1 = x;
          yh = y;
        }

        if (x == 1.0 || x == 0.0)
        {
          break;
        }

        // Compute the derivative of the function at this point
        d = ((a - 1.0) * Math.Log(x)) + ((b - 1.0) * Math.Log(1.0 - x)) + lgm;

        if (d < Constants.MinimumLogResults)
        {
          goto done;
        }

        if (d > Constants.MaximumLogResults)
        {
          break;
        }

        //! #8584 check that the evaluation of (y-y0)/Exp(d) will not overflow
        if (y != y0 && (Math.Log(Math.Abs(y - y0)) - d) > Constants.MaximumLogResults)
        {
          break;
        }

        d = Math.Exp(d);
        // Compute the step to the next approximation of x
        d = (y - y0) / d;
        double xt = x - d;
        if (xt <= x0)
        {
          y = (x - x0) / (x1 - x0);
          xt = x0 + (0.5 * y * (x - x0));
          if (xt <= 0.0)
          {
            break;
          }
        }

        if (xt >= x1)
        {
          y = (x1 - x) / (x1 - x0);
          xt = x1 - (0.5 * y * (x1 - x));
          if (xt >= 1.0)
          {
            break;
          }
        }

        x = xt;
        if (CheckRatio(d, x, newtonTolerance))
        {
          goto done;
        }
      }

      //! Did not converge

      dithresh = 2 * newtonTolerance;
      goto ihalve;

    //!+ done:
    done:

      if (isReversed)
      {
        if (x <= Constants.MachinePrecision)
        {
          x = 1.0 - Constants.MachinePrecision;
        }
        else
        {
          x = 1.0 - x;
        }
      }

      return x;
    }
    [MethodImpl(MethodImplOptions.AggressiveInlining)] // Though one would hope this inlines anyway
    private static bool CheckRatio(double numerator, double denominator, double limit)
    {
      return Math.Abs(numerator) < Math.Abs(denominator * limit);
    }
    [MethodImpl(MethodImplOptions.AggressiveInlining)] // Though one would hope this inlines anyway
    private static bool CheckRatioPositiveDenom(double numerator, double positiveDenominator, double limit)
    {
      return Math.Abs(numerator) < positiveDenominator * limit;
    }


    /// <summary>
    ///   The Log Gamma function.
    /// </summary>
    /// <param name="x">The parameter to calculate the Log Gamma of.</param>
    /// <returns>The natural logarithm of Gamma(<paramref name="x" />).</returns>
    /// <exception cref="ArithmeticException">Thrown if an overflow occurs.</exception>
    /// <remarks>
    ///   <para>This is the Cephes <c>lgam</c> method.</para>
    ///   <para>Note that Log Gamma(x + 1) ≈ Log(x!).</para>
    ///   <para>The true result is complex for x &lt; 0.</para>
    /// </remarks>
    internal static double LogGamma(double x)
    {
      const double cMaxLogGamma = 2.556348e305;
      const double cLogSqrt2Pi = 0.91893853320467274178; // log(sqrt(2 * pi))

      if (x < -34)
      {
        return LogGammaVerySmallX(x);
      }

      if (x < 13)
      {
        return LogGammaSmallX(x);
      }

      if (cMaxLogGamma < x)
      {
        throw new ArithmeticException("Overflow calculating the log gamma");
      }

      double q = ((x - 0.5) * Math.Log(x)) - x + cLogSqrt2Pi;

      if (1.0e8 < x)
      {
        return q;
      }

      double p = 1 / (x * x);

      if (1000 <= x)
      {
        q += ((((7.9365079365079365079365e-4 * p)
                - 2.7777777777777777777778e-3) * p)
              + 0.0833333333333333333333) / x;
      }
      else
      {
        q += PolyEval_1(p, cA) / x;
      }

      return q;
    }

    /// <summary>
    ///   Calculates the unscaled inverse Normal.
    /// </summary>
    /// <param name="u">The probability.</param>
    /// <returns>The inverse of Normal(0, 1) for the supplied probability.</returns>
    /// <remarks>
    ///   <para>Not actually from Cephes, but this seems the obvious place to put it.</para>
    ///   <para>This function is also used for calculating Poisson.</para>
    /// </remarks>
    internal static double MorosInvNormal(double u)
    {
      double r;

      const double a1 = 2.50662823884; // sqrt 2 pi?
      const double a2 = -18.61500062529;
      const double a3 = 41.39119773534;
      const double a4 = -25.44106049637;
      const double b1 = -8.4735109309;
      const double b2 = 23.08336743743;
      const double b3 = -21.06224101826;
      const double b4 = 3.13082909833;
      const double c1 = 0.337475482272615;
      const double c2 = 0.976169019091719;
      const double c3 = 0.160797971491821;
      const double c4 = 2.76438810333863E-02;
      const double c5 = 3.8405729373609E-03;
      const double c6 = 3.951896511919E-04;
      const double c7 = 3.21767881768E-05;
      const double c8 = 2.888167364E-07;
      const double c9 = 3.960315187E-07;

      double x = u - 0.5;

      if (Math.Abs(x) < 0.42)
      {
        r = x * x;
        r = x * ((((((a4 * r) + a3) * r) + a2) * r) + a1) / ((((((((b4 * r) + b3) * r) + b2) * r) + b1) * r) + 1);
      }
      else
      {
        if (0 < x)
        {
          r = Math.Log(-Math.Log(1 - u));
        }
        else
        {
          r = Math.Log(-Math.Log(u));
        }

        r = c1 + (r * (c2 + (r * (c3 + (r * (c4 + (r * (c5 + (r * (c6 + (r * (c7 + (r * (c8 + (r * c9)))))))))))))));

        if (x <= 0)
        {
          r = -r;
        }
      }

      return r;
    }

    /// <summary>
    ///   Calculates the complement of the incomplete Gamma function (upper).
    /// </summary>
    /// <param name="a">The a parameter.</param>
    /// <param name="x">The x parameter.</param>
    /// <returns>
    ///   The complement of the incomplete Gamma function.
    /// </returns>
    /// <remarks>
    ///   This is the Cephes <c>igamc</c> method.
    /// </remarks>
    private static double ComplementedIncompleteGamma(double a, double x)
    {
      if (Double.IsNaN(a) || Double.IsNaN(x))
      {
        throw new ArithmeticException("Floating point overflow");
      }

      if (x <= 0 || a <= 0)
      {
        return 1;
      }

      if (x < 1 || x < a)
      {
        return 1 - IncompleteGamma(a, x);
      }

      double ax = (a * Math.Log(x)) - x - LogGamma(a);

      if (ax < -Constants.MaximumLogResults)
      {
        return 0;
      }

      ax = Math.Exp(ax);
      if (double.IsInfinity(ax))
      {
        // we can't return a value in [0, 1] if ax is infinity
        throw new ArithmeticException("Floating point overflow");
      }

      //! continued fraction

      double y = 1 - a;
      double z = x + y + 1;
      double c = 0;
      double pkm2 = 1;
      double qkm2 = x;
      double pkm1 = x + 1;
      double qkm1 = z * x;
      double ans = pkm1 / qkm1;

      do
      {
        c += 1;
        y += 1;
        z += 2;
        double yc = y * c;

        double pk = (pkm1 * z) - (pkm2 * yc);
        if (double.IsInfinity(pk))
        {
          throw new ArithmeticException("Floating point overflow");
        }

        double qk = (qkm1 * z) - (qkm2 * yc);
        if (double.IsInfinity(qk))
        {
          throw new ArithmeticException("Floating point overflow");
        }

        if (qk != 0)
        {
          double r = pk / qk;
          if (CheckRatio(ans - r, r, Constants.MachinePrecision))
          {
            break;
          }
          ans = r;
        }

        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if (cBig < Math.Abs(pk))
        {
          pkm2 *= cBigInv;
          pkm1 *= cBigInv;
          qkm2 *= cBigInv;
          qkm1 *= cBigInv;
        }
      } while (true);

      double result = ans * ax;

      if (result < 0 || 1 < result)
      {
        string message = string.Format(CultureInfo.InvariantCulture, "Invalid return value {0:R} in ComplementedIncompleteGamma", result);
        throw new ArithmeticException(message);
      }

      return result;
    }

    /// <summary>
    ///   Calculates the gamma function at the specified value of x.
    /// </summary>
    /// <exception cref="UserArgumentException">Thrown when a user argument error condition occurs.</exception>
    /// <param name="x">The x parameter.</param>
    /// <returns>
    ///   The result of the gamma function at the specified value of x.
    /// </returns>
    /// <remarks>
    ///   This is the <c>gam</c> function in Cephes.
    /// </remarks>
    internal static double Gamma(double x)
    {
      double xOnEntry = x;

      int sgngam = 1;
      if (double.IsNaN(x))
      {
        return x;
      }

      if (double.IsPositiveInfinity(x))
      {
        return x;
      }

      if (double.IsNegativeInfinity(x))
      {
        return double.NaN;
      }

      double p;
      double q = Math.Abs(x);

      double z;

      if (q > 33.0)
      {
        if (x < 0.0)
        {
          p = Math.Truncate(q);
          if (p == q)
          {
            //Debugger.Launch();
            throw new UserArgumentException($"Overflow(1) in Gamma function: x={x}");
          }

          int i = (int)Math.Truncate(p);
          if ((i & 0x1) == 0)
          {
            sgngam = -1;
          }

          z = q - p;
          if (z > 0.5)
          {
            p = p + 1.0;
            z = q - p;
          }

          z = q * Math.Sin(Math.PI * z);
          if (z == 0.0)
          {
            //Debugger.Launch();
            throw new UserArgumentException($"Overflow(2) in Gamma function: x={x}");
          }

          z = Math.Abs(z);
          z = Math.PI / (z * GammaUsingStirlingsFormula(q));
        }
        else
        {
          z = GammaUsingStirlingsFormula(x);
        }

        return sgngam * z;
      }

      z = 1.0;
      while (x >= 3.0)
      {
        x = x - 1.0;
        z = z * x;
      }

      while (x < 0.0)
      {
        if (x > -1E-9)
        {
          goto small;
        }

        z = z / x;
        x = x + 1.0;
      }

      while (x < 2.0)
      {
        if (x < 1e-9)
        {
          goto small;
        }

        z = z / x;
        x = x + 1.0;
      }

      if (x == 2.0)
      {
        return z;
      }

      x = x - 2.0;
      p = PolyEval(x, cPgam);
      q = PolyEval(x, cQgam);
      return z * p / q;

    //!+ small:
    small:
      if (x == 0.0)
      {
        //Debugger.Launch();
        throw new UserArgumentException($"Overflow(3) in Gamma function: x={xOnEntry} on entry");
      }

      return z / ((1.0 + (0.5772156649015329 * x)) * x);
    }

    /// <summary>
    ///   Calculates the gamma function using Stirling's formula.
    /// </summary>
    /// <param name="x">The x parameter. The polynomial STIR is valid for 33 &lt;= x &lt;= 172.</param>
    /// <returns>
    ///   The value of the gamma function at the specified value of x.
    /// </returns>
    /// <remarks>
    ///   This is the <c>stirf</c> function in Cephes.
    /// </remarks>
    private static double GammaUsingStirlingsFormula(double x)
    {
      const double cMaxStirling = 143.01608;
      const double cSqrt2Pi = 2.50662827463100050242E0; // sqrt(2 *pi)

      double w = 1.0 / x;
      w = 1.0 + (w * PolyEval(w, cStirling));
      double yInv = Math.Exp(-x);

      if (x > cMaxStirling)
      {
        // Avoid overflow in pow()
        double v = Math.Pow(x, (0.5 * x) - 0.25);
        return cSqrt2Pi * w * v * (v * yInv);
      }
      else
      {
        return cSqrt2Pi * w * Math.Pow(x, x - 0.5) * yInv;
      }
    }

    /// <summary>
    ///   Calculates the incomplete beta integral using continued fraction expansion (#1).
    /// </summary>
    /// <param name="a">The alpha parameter.</param>
    /// <param name="b">The beta parameter.</param>
    /// <param name="x">The x parameter.</param>
    /// <param name="tolerance">The required accuracy for iterative convergence. 0=> use default (3 epsilon)</param>
    /// <returns>
    ///   The value of the incomplete beta integral with the supplied parameters.
    /// </returns>
    /// <remarks>
    ///   This is the <c>incbcf</c> function in Cephes.
    /// </remarks>
    private static double IncompleteBetaIntegral1(double a, double b, double x, Tolerance tolerance)
    {
      double k1 = a;
      double k2 = a + b;
      double k3 = a;
      double k4 = a + 1.0;
      double k5 = 1.0;
      double k6 = b - 1.0;
      double k7 = k4;
      double k8 = a + 2.0;

      double pkm2 = 0.0;
      double qkm2 = 1.0;
      double pkm1 = 1.0;
      double qkm1 = 1.0;
      double ans = 1.0;
      double r = 1.0;
      int n = 0;
      double thresh = tolerance.Value == 0 ? 3.0 * Constants.MachinePrecision : tolerance.Value;

      do
      {
        double xk = -(x * k1 * k2) / (k3 * k4);
        double pk = pkm1 + (pkm2 * xk);
        double qk = qkm1 + (qkm2 * xk);
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        xk = (x * k5 * k6) / (k7 * k8);
        pk = pkm1 + (pkm2 * xk);
        qk = qkm1 + (qkm2 * xk);
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        if (qk != 0)
        {
          r = pk / qk;
        }

        if (r != 0)
        {
          // was comparing |(ans-r)/r|<thresh, but since r is positive, we can use a multiply instead of a divide
          if (CheckRatioPositiveDenom(ans - r, r, thresh))
          {
            ans = r;
            break;
          }
        }

        ans = r;
        k1 = k1 + 1.0;
        k2 = k2 + 1.0;
        k3 = k3 + 2.0;
        k4 = k4 + 2.0;
        k5 = k5 + 1.0;
        k6 = k6 - 1.0;
        k7 = k7 + 2.0;
        k8 = k8 + 2.0;

        if ((Math.Abs(qk) + Math.Abs(pk)) > cBig)
        {
          pkm2 *= cBigInv;
          pkm1 *= cBigInv;
          qkm2 *= cBigInv;
          qkm1 *= cBigInv;
        }

        else if ((Math.Abs(qk) < cBigInv) || (Math.Abs(pk) < cBigInv))
        {
          pkm2 *= cBig;
          pkm1 *= cBig;
          qkm2 *= cBig;
          qkm1 *= cBig;
        }

        n++;
      } while (n <= 300);

      return ans;
    }

    /// <summary>
    ///   Calculates the incomplete beta integral using continued fraction expansion (#2).
    /// </summary>
    /// <param name="a">The alpha parameter.</param>
    /// <param name="b">The beta parameter.</param>
    /// <param name="x">The x parameter.</param>
    /// <param name="tolerance">The required accuracy for terminating the iteration. 0=> default of 3*epsilon</param>
    /// <returns>
    ///   The value of the incomplete beta integral with the supplied parameters.
    /// </returns>
    /// <remarks>
    ///   This is the <c>incbd</c> function in Cephes.
    /// </remarks>
    private static double IncompleteBetaIntegral2(double a, double b, double x, Tolerance tolerance)
    {
      double k1 = a;
      double k2 = b - 1.0;
      double k3 = a;
      double k4 = a + 1.0;
      double k5 = 1.0;
      double k6 = a + b;
      double k7 = a + 1.0;
      double k8 = a + 2.0;

      double pkm2 = 0.0;
      double qkm2 = 1.0;
      double pkm1 = 1.0;
      double qkm1 = 1.0;
      double z = x / (1.0 - x);

      double ans = 1.0;
      double r = 1.0;
      int n = 0;
      double thresh = tolerance.Value == 0 ? 3.0 * Constants.MachinePrecision : tolerance.Value;
      do
      {
        double xk = -(z * k1 * k2) / (k3 * k4);
        double pk = pkm1 + (pkm2 * xk);
        double qk = qkm1 + (qkm2 * xk);
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        xk = (z * k5 * k6) / (k7 * k8);
        pk = pkm1 + (pkm2 * xk);
        qk = qkm1 + (qkm2 * xk);
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        if (qk != 0)
        {
          r = pk / qk;
        }


        if (r != 0)
        {

          if (CheckRatio(ans - r, r, thresh))
          {
            ans = r;
            return ans;
          }

          ans = r;
          if (double.IsNaN(ans))
          {
            // it's never going to become non-NaN again
            return ans;
          }
        }


        k1 = k1 + 1.0;
        k2 = k2 - 1.0;
        k3 = k3 + 2.0;
        k4 = k4 + 2.0;
        k5 = k5 + 1.0;
        k6 = k6 + 1.0;
        k7 = k7 + 2.0;
        k8 = k8 + 2.0;

        if (Math.Abs(qk) + Math.Abs(pk) > cBig)
        {
          pkm2 = pkm2 * cBigInv;
          pkm1 = pkm1 * cBigInv;
          qkm2 = qkm2 * cBigInv;
          qkm1 = qkm1 * cBigInv;
        }

        else if (Math.Abs(qk) < cBigInv || Math.Abs(pk) < cBigInv)
        {
          pkm2 = pkm2 * cBig;
          pkm1 = pkm1 * cBig;
          qkm2 = qkm2 * cBig;
          qkm1 = qkm1 * cBig;
        }

        n++;
      } while (n < 300);

      return ans;
    }

    public static double lastAPlusB, lastA, lastB;

    /// <summary>
    ///   Calculates the power series for the incomplete beta integral.
    /// </summary>
    /// <param name="a">The alpha parameter.</param>
    /// <param name="b">The beta parameter.</param>
    /// <param name="tolerance">Allowable relative error or 0 to use default</param>
    /// <param name="x">The x parameter.</param>
    /// <returns>
    ///   The result of the power series.
    /// </returns>
    /// <remarks>
    ///   <para>Use when b*x is small and x not too close to 1.</para>
    ///   <para>This is the <c>pseries</c> function in Cephes.</para>
    /// </remarks>
    private static double IncompleteBetaIntegralPowerSeries(double a, double b, Tolerance tolerance, double x)
    {
      double toleranceVal = tolerance.Value == 0 ? Constants.MachinePrecision : tolerance.Value;
      double s = 0.0;
      double ai = 1.0 / a;
      double t1;
      {

        double u2 = (1.0 - b) * x;
        double v = u2 / (a + 1.0);
        t1 = v;
        double t2 = u2;
        double n = 2.0;
        double z = toleranceVal * ai;
        double one = 1.0;
        s = IterateForS(a, b, x, v, t2, n, s, z, one);
      }
      s = s + t1;
      s = s + ai;

      double u = a * Math.Log(x);
      double t;
      double aplusb = a + b;
      if (aplusb < cMaxGam && Math.Abs(u) < Constants.MaximumLogResults)
      {
        lastAPlusB = aplusb;
        lastA = a;
        lastB = b;
        double gammaAPlusB = Gamma(aplusb);
        double gammaA = Gamma(a);
        double gammaB = Gamma(b);

        t = gammaAPlusB / (gammaA * gammaB);

        //! If this is zero and gammaAPlusB != 0, see if we can do better with logs
        //! Without this correction, you can get an infinite loop
        if (t == 0 && gammaAPlusB != 0)
        {
          t = Math.Exp(Math.Log(gammaAPlusB) - Math.Log(gammaA) - Math.Log(gammaB));
        }

        s = s * t * Math.Exp(u); //pow(x,a);
      }
      else
      {
        //t = LogGamma(aplusb) - LogGamma(a) - LogGamma(b) + u + Math.Log(s);
        if (aplusb == a)
        {
          t = KahanSum(-LogGamma(b), u, Math.Log(s));
        }
        else if (aplusb == b)
        {
          t = KahanSum(-LogGamma(a), u, Math.Log(s));
        }
        else
        {
          t = KahanSum(LogGamma(aplusb), -LogGamma(a), -LogGamma(b), u, Math.Log(s));
        }

        s = Math.Exp(t);
      }

      return s;
    }
    /// <summary>
    /// Having this as a separate function allows for better register allocation
    /// </summary>
    private static double IterateForS(double a, double b, double x, double v, double t, double n, double s, double z, double one)
    {
      while (Math.Abs(v) > z)
      {
        t = t * (n - b) * (x / n);
        v = t / (a + n);
        s = s + v;
        n = n + one;
      }
      return s;
    }

    public static double lastAForPowerSeries;

    /// <summary>
    ///   Calculates the incomplete beta integral.
    /// </summary>
    /// <exception cref="UserArgumentException">Thrown when an error condition occurs.</exception>
    /// <param name="aa">The alpha parameter.</param>
    /// <param name="bb">The beta parameter.</param>
    /// <param name="tolerance">The allowed relative tolerance or 0 for default</param>
    /// <param name="xx">The x parameter.</param>
    /// <returns>
    ///   Returns incomplete beta integral of the arguments, evaluated from zero to x.
    /// </returns>
    /// <remarks>
    ///   This is <c>incbet</c> in the Delphi code.
    /// </remarks>
    internal static double IncompleteBetaIntegral(double aa, double bb, Tolerance tolerance, double xx)
    {
      const double cCephesPeakRelativeError = 8.7e-10;

      if (aa <= 0.0 || bb <= 0.0)
      {
        //Debugger.Launch();
        throw new UserArgumentException("Both parameters for Beta function must be positive");
      }

      if (xx <= 0.0 || xx >= 1.0)
      {
        if (xx == 0.0)
        {
          return 0;
        }

        if (xx == 1.0)
        {
          return 1;
        }

        throw new UserArgumentException("Argument for Beta function must be non-negative");
      }

      bool isReversed;

      double a;
      double b;
      double xc;
      double x;

      double t;

      if (bb * xx <= 1.0 && xx <= 0.95)
      {
        t = IncompleteBetaIntegralPowerSeries(aa, bb, tolerance, xx);

        //! This goto prevented a and b from being initialised, which the compiler correctly reported.
        isReversed = false;
        a = aa;
        b = bb;

        goto done;
      }

      double w = 1.0 - xx;

      // Reverse a and b if x is greater than the mean
      if (xx > aa / (aa + bb))
      {
        isReversed = true;
        a = bb;
        b = aa;
        xc = xx;
        x = w;
      }
      else
      {
        isReversed = false;
        a = aa;
        b = bb;
        xc = w;
        x = xx;
      }

      if (isReversed && b * x <= 1.0 && x <= 0.95)
      {
        lastAForPowerSeries = a;
        t = IncompleteBetaIntegralPowerSeries(a, b, tolerance, x);
        goto done;
      }

      // Choose expansion for better convergence
      double y = (x * (a + b - 2.0)) - (a - 1.0);
      if (y < 0.0)
      {
        w = IncompleteBetaIntegral1(a, b, x, tolerance);
      }
      else
      {
        w = IncompleteBetaIntegral2(a, b, x, tolerance) / xc;
      }

      // Multiply w by the factor
      //     a      b   _             _     _
      //    x  (1-x)   | (a+b) / ( a | (a) | (b) )

      y = a * Math.Log(x);
      t = b * Math.Log(xc);
      double aplusb = a + b;
      if (aplusb < cMaxGam && Math.Abs(y) < Constants.MaximumLogResults && Math.Abs(t) < Constants.MaximumLogResults)
      {
        t = Math.Exp(t) * (Math.Exp(y) / a) * w * (Gamma(aplusb) / (Gamma(a) * Gamma(b)));
        goto done;
      }

      // Resort to logarithms
      double logGammaAPlusB = LogGamma(a + b);
      double logGammaA = LogGamma(a);
      double logGammaB = LogGamma(b);
      double logWOverA = Math.Log(w / a);

      y = y + t + logGammaAPlusB - logGammaA - logGammaB + logWOverA;
      //x  y = y + t + LogGamma(a + b) - LogGamma(a) - LogGamma(b);
      //x  y = y + Math.Log(w / a);

      if (y < Constants.MinimumLogResults)
      {
        t = 0.0;
      }
      else
      {
        t = Math.Exp(y);
      }

    //!+ done:
    done:

      if (isReversed)
      {
        if (t <= Constants.MachinePrecision)
        {
          t = 1.0 - Constants.MachinePrecision;
        }
        else
        {
          t = 1.0 - t;
        }
      }

      double result = t;

      // Check for a result that is greater than 1.0 or less than 0.0 (invalid).
      if (result < 0.0 || result > 1.0)
      {
        if (result > 1.0 + cRootMachinePrecision || result < -cRootMachinePrecision)
        {
          //Check for a case where the result should be close to zero.
          if (b > 1.0)
          {
            // This may be due to a very small xx and a very large bb, in which case the
            // answer should be very close to zero.  We can check for this case by
            // computing an upper bound on the result in terms of xx and aa.
            // First, compute an upper bound on the numerator of the incomplete beta
            // function for this case:
            double upperBound = Math.Pow(xx, aa) / aa;
            Debug.Assert(0 <= upperBound, "Upper bound should be non-negative");

            //Divide the upper bound by a Beta function.  This
            //represents the denominator in the expression for the incomplete Beta
            //function:
            upperBound = upperBound * Math.Exp(LogGamma(aa + bb) - LogGamma(aa) - LogGamma(bb));

            //Cephes quotes a peak relative error of 8.7e-10.  Thus, if the upper bound
            //is below twice this then return half of the the upper bound as the answer.
            if (upperBound < 2.0 * cCephesPeakRelativeError)
            {
              //We know that the result is in [0, UpperBound], so return the middle of
              //this range.
              return upperBound / 2.0;
            }
          }

          //Check for a case where the result should be close to one.  This uses essentially
          //the same logic as the above case, and follows from the symmetry relationship
          //Incbet(x)(a, b) = 1 - Incbet(1-x)(b, a).
          if (a > 1.0)
          {
            // This may be due to a xx close to one and a very large aa, in which case the
            // answer should be very close to one.  We can check for this case by
            // computing a lower bound on the result in terms of xx and bb, via
            // computing an upper bound on Incbet(1-x)(b, a).

            double upperBound = Math.Pow(1.0 - xx, bb) / bb;
            Debug.Assert(0 <= upperBound, "Upper bound should be non-negative");
            upperBound = upperBound * Math.Exp(LogGamma(aa + bb) - LogGamma(aa) - LogGamma(bb));

            //Cephes quotes a peak relative error of 8.7e-10.  Thus, if the upper bound
            //is below twice this then return (1.0 - half of the the upper bound) as the answer.
            if (upperBound < 2.0 * cCephesPeakRelativeError)
            {
              //We know that the result is in [1.0 - UpperBound, 1.0], so return the middle of
              //this range.
              return 1.0 - (upperBound / 2.0);
            }
          }

          //          if (!TAlaskaConfiguration.SuppressIncBetValidation)
          //          {
          throw new UserArgumentException("Cannot compute an answer for the incremental beta function with the given arguments (x={0}, a={1}, b={2}). Proposed result was {3}.", xx,
            aa, bb, result);
          //          }
        }

        if (result < 0.0)
        {
          // Result is only just outside of the range [0, 1]
          Debug.Assert(result >= -cRootMachinePrecision, "Result should be only just negative");

          return 0;
        }

        Debug.Assert(1 < result && result <= 1 + cRootMachinePrecision, "We're here because the result is (slightly) out of bounds");

        return 1;
      }

      return result;
    }

    /*							Stdtr.c
     *
     *	Student's t distribution
     *
     *
     *
     * SYNOPSIS:
     *
     * double t, Stdtr();
     * short k;
     *
     * y = Stdtr( k, t );
     *
     *
     * DESCRIPTION:
     *
     * Computes the integral from minus infinity to t of the Student
     * t distribution with integer k > 0 degrees of freedom:
     *
     *                                      t
     *                                      -
     *                                     | |
     *              -                      |         2   -(k+1)/2
     *             | ( (k+1)/2 )           |  (     x   )
     *       ----------------------        |  ( 1 + --- )        dx
     *                     -               |  (      k  )
     *       sqrt( k pi ) | ( k/2 )        |
     *                                   | |
     *                                    -
     *                                   -inf.
     *
     * Relation to incomplete beta integral:
     *
     *        1 - Stdtr(k,t) = 0.5 * incbet( k/2, 1/2, z )
     * where
     *        z = k/(k + t**2).
     *
     * For t < -2, this is the method of computation.  For higher t,
     * a direct method is derived from integration by parts.
     * Since the function is symmetric about t=0, the area under the
     * right tail of the density is found by calling the function
     * with -t instead of t.
     *
     * ACCURACY:
     *
     * Tested at random 1 <= k <= 25.  The "domain" refers to t.
     *                      Relative error:
     * arithmetic   domain     # trials      peak         rms
     *    IEEE     -100,-2      50000       5.9e-15     1.4e-15
     *    IEEE     -2,100      500000       2.7e-15     4.9e-17
    */


    /// <summary>
    ///   Calculates the Log Gamma function when X is very small.
    /// </summary>
    /// <param name="x">The parameter to calculate the Log Gamma of.</param>
    /// <returns>The natural logarithm of Gamma(<paramref name="x" />).</returns>
    private static double LogGammaSmallX(double x)
    {
      Debug.Assert(-34 <= x);
      Debug.Assert(x < 13);

      double z = 1;
      double p = 0;
      double u = x;

      while (3 <= u)
      {
        p--;
        u = x + p;
        z *= u;
      }

      while (u < 2)
      {
        if (u == 0)
        {
          throw new ArithmeticException("Overflow calculating log gamma for a small value of x");
        }

        z /= u;
        p++;
        u = x + p;
      }

      if (z < 0)
      {
        z = -z;
      }

      if (u == 2)
      {
        return Math.Log(z);
      }

      p -= 2;
      x += p;
      p = x * PolyEval(x, cB) / PolyEval_1(x, cC);

      double result = Math.Log(z) + p;

      return result;
    }

    /// <summary>
    ///   Calculates the Log Gamma function when X is very small.
    /// </summary>
    /// <param name="x">The parameter to calculate the Log Gamma of.</param>
    /// <returns>The natural logarithm of Gamma(<paramref name="x" />).</returns>
    private static double LogGammaVerySmallX(double x)
    {
      Debug.Assert(x < -34);

      const double cLogPi = 1.14472988584940017414;

      double q = -x;
      double w = LogGamma(q);
      double p = Math.Floor(q);
      if (p == q)
      {
        throw new ArithmeticException("Overflow calculating log gamma for a very small value of x");
      }

      double z = q - p;

      if (0.5 < z)
      {
        p++;
        z = p - q;
      }

      z = q * Math.Sin(Math.PI * z);

      if (z == 0)
      {
        throw new ArithmeticException("Overflow calculating log gamma for a very small value of x (#2)");
      }

      z = cLogPi - Math.Log(z) - w;

      return z;
    }

    /// <summary>
    ///   <para>Evaluates the polynomial described by <paramref name="coefficients" /> at point <paramref name="x" />.</para>
    ///   <para>For example, PolyEval(x, [c₀, c₁, c₂]) = c₀x² + c₁x + c₂.</para>
    /// </summary>
    /// <param name="x">The value at which to evaluate the polynomial.</param>
    /// <param name="coefficients">The coefficients of the polynomial in descending order.</param>
    /// <returns>The value of the polynomial described by <paramref name="coefficients" /> at point <paramref name="x" />.</returns>
    /// <remarks>Unlike Delphi routine <c>Poly</c>, this assumes that coefficients array is in descending order.</remarks>
    internal static double PolyEval(double x, double[] coefficients)
    {
      double result = coefficients[0];

      for (int i = 1; i < coefficients.Length; i++)
      {
        result = (result * x) + coefficients[i];
      }

      return result;
    }

    /// <inheritdoc cref="PolyEval" />
    /// <summary>
    ///   <para>Evaluates the polynomial described by <paramref name="coefficients" /> at point <paramref name="x" />.</para>
    ///   <para>For example, PolyEval_1(x, [c₀, c₁, c₂]) = x³ + c₀x² + c₁x + c₂.</para>
    /// </summary>
    private static double PolyEval_1(double x, double[] coefficients)
    {
      double result = x;

      for (int i = 0; i < coefficients.Length; i++)
      {
        result = (result * x) + coefficients[i];
      }

      return result;
    }

    private static readonly double[] cPerf =
    {
      2.46196981473530512524E-10,
      5.64189564831068821977E-1,
      7.46321056442269912687E0,
      4.86371970985681366614E1,
      1.96520832956077098242E2,
      5.26445194995477358631E2,
      9.34528527171957607540E2,
      1.02755188689515710272E3,
      5.57535335369399327526E2
    };

    private static readonly double[] cQerf =
    {
// 1.00000000000000000000E0,
      1.32281951154744992508E1,
      8.67072140885989742329E1,
      3.54937778887819891062E2,
      9.75708501743205489753E2,
      1.82390916687909736289E3,
      2.24633760818710981792E3,
      1.65666309194161350182E3,
      5.57535340817727675546E2
    };

    private static readonly double[] cRerf =
    {
      5.64189583547755073984E-1,
      1.27536670759978104416E0,
      5.01905042251180477414E0,
      6.16021097993053585195E0,
      7.40974269950448939160E0,
      2.97886665372100240670E0
    };

    private static readonly double[] cSerf =
    {
// 1.00000000000000000000E0,
      2.26052863220117276590E0,
      9.39603524938001434673E0,
      1.20489539808096656605E1,
      1.70814450747565897222E1,
      9.60896809063285878198E0,
      3.36907645100081516050E0
    };


    private static readonly double[] cTerf =
    {
      9.60497373987051638749E0,
      9.00260197203842689217E1,
      2.23200534594684319226E3,
      7.00332514112805075473E3,
      5.55923013010394962768E4
    };

    private static readonly double[] cUerf =
    {
      // 1.00000000000000000000E0,
      3.35617141647503099647E1,
      5.21357949780152679795E2,
      4.59432382970980127987E3,
      2.26290000613890934246E4,
      4.92673942608635921086E4
    };


    public static double Erf(double x)

    {
      // THREAD-SAFE

      if (Math.Abs(x) > 1.0)
      {
        return 1.0 - Erfc(x);
      }

      double z = x * x;
      double y = x * PolyEval(z, cTerf) / PolyEval_1(z, cUerf);
      return y;
    }


    private const double cMaxM = 128.0;
    private const double cMinv = 0.0078125;
    public static readonly double cMaxlog = Math.Log(double.MaxValue);

    public static double Expx2(double x, int sign)
    {
      // THREAD-SAFE
      x = Math.Abs(x);
      if (sign < 0)
      {
        x = -x;
      }

      /*Represent x as an exact multiple of M plus a residual.
        M is a power of 2 chosen so that exp(m * m) does not overflow
        or underflow and so that | x - m | is small.  */
      double m = cMinv * Math.Truncate(cMaxM * x + 0.5);
      double f = x - m;

      // x^2 = m^2 + 2mf + f^2 */
      double u = m * m;
      double u1 = 2 * m * f + f * f;

      if (sign < 0)
      {
        u = -u;
        u1 = -u1;
      }

      if ((u + u1) > cMaxlog)
      {
        return double.PositiveInfinity;
      }

      // u is exact, u1 is small.
      return Math.Exp(u) * Math.Exp(u1);
    }

    public static double Erfc(double a)

    {
      double x = a < 0.0 ? -a : a;

      if (x < 1.0)
      {
        return 1.0 - Erf(a); // THREAD-SAFE
      }

      double z = -a * a;

      if ((z < -cMaxlog))
      {
        if (a < 0)
        {
          return 2.0;
        }

        return 0.0;
      }

      // Compute z = exp(z).  */
      z = Expx2(a, -1);
      double p;
      double q;
      if (x < 8.0)
      {
        p = PolyEval(x, cPerf); // THREAD-SAFE
        q = PolyEval_1(x, cQerf); // THREAD-SAFE
      }
      else
      {
        p = PolyEval(x, cRerf); // THREAD-SAFE
        q = PolyEval_1(x, cSerf); // THREAD-SAFE
      }

      double y = (z * p) / q;

      if (a < 0)
      {
        y = 2.0 - y;
      }

      return y;
    }


    /* Exponentially scaled Erfc function
       exp(x^2) Erfc(x)
       valid for x > 1.
       Use with Ndtr and Expx2.  */

    static double KahanSum(params double[] args)
    {
      double sum = 0.0;
      double c = 0.0;
      foreach (double d in args)
      {
        double y = d - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
      }

      return sum;
    }
  }
}