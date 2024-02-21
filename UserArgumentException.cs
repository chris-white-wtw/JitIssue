using System.Globalization;

namespace JitIssue
{
  /// <summary>
  ///   This class represents an exception arising from bad arguments supplied to a maths formula, which
  ///   should be reported back to the user (i.e. the description should make sense to them).
  /// </summary>
  [Serializable]
  public class UserArgumentException : ArithmeticException
  {
    #region Constructors

    /// <summary>
    ///   Initialises a new instance of the <see cref="UserArgumentException" /> class.
    /// </summary>
    public UserArgumentException()
    {
    }

    /// <summary>
    ///   Initialises a new instance of the <see cref="UserArgumentException" /> class with a specified error message.
    /// </summary>
    /// <param name="message">The error message that explains the reason for the exception.</param>
    public UserArgumentException(string message)
      : base(message)
    {
    }

    /// <summary>
    ///   Initialises a new instance of the <see cref="UserArgumentException" /> class with a specified error message and a
    ///   reference to the inner exception that is the cause of this exception.
    /// </summary>
    /// <param name="message">The error message that explains the reason for the exception.</param>
    /// <param name="innerException">
    ///   The exception that is the cause of the current exception. If the
    ///   <paramref name="innerException" /> parameter is not a null reference, the current exception is raised in a catch
    ///   block that handles the inner exception.
    /// </param>
    public UserArgumentException(string message, Exception innerException)
      : base(message, innerException)
    {
    }

    /// <summary>
    ///   Initialises a new instance of the <see cref="UserArgumentException" /> class with a format string and some details to
    ///   populate the format string with.
    /// </summary>
    /// <param name="messageFormat">
    ///   The formatted error message. This should include numbered placeholders as used by
    ///   <see cref="String.Format(string, object[])" />
    /// </param>
    /// <param name="messageContents">Objects supplying supplemental information to be inserted into the error message.</param>
    public UserArgumentException(string messageFormat, params object[] messageContents)
      : base(string.Format(CultureInfo.InvariantCulture, messageFormat, messageContents))
    {
    }

    #endregion Constructors
  }
}