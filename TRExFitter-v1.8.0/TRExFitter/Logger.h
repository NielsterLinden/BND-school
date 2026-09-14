/**
 * @file Logger.h
 * @brief Logging class
 *
 */

#pragma once

#include <chrono>
#include <cstring>
#include <ctime>
#include <iostream>

#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#define LOG(x) Logger::get()(LoggingLevel::x, __FILENAME__, __LINE__, __FUNCTION__)

/**
 * @brief Enum storing different logging levels
 *
 */
enum class LoggingLevel {
  ERROR = 0,
  WARNING = 1,
  INFO = 2,
  DEBUG = 3,
  VERBOSE = 4
};

/**
 * @brief Singletop class for logging
 *
 */
class Logger {
public:

  /**
   * @brief Deleted copy constructor
   *
   */
  Logger(const Logger&) = delete;

  /**
   * @brief Deleted assignment operator
   *
   */
  void operator=(const Logger&) = delete;

  /**
   * @brief Returns the class. Created on first call, then persistent
   *
   * @return Logger&
   */
  static Logger& get() {
    static Logger logger;
    return logger;
  }

  /**
   * @brief Set the Log Level object
   *
   * @param level
   */
  void setLogLevel(const LoggingLevel& level) {
    m_logLevel = level;
  }

  /**
   * @brief Set whether to include the lines in the logging
   *
   * @param flag
   */
  void setIncludeLines(const bool flag) {
    m_includeLines = flag;
  }

  /**
   * @brief Set whether to include function in the log
   *
   * @param flag
   */
  void setIncludeFunction(const bool flag) {
    m_includeFunction = flag;
  }

  /**
   * @brief Set whether to include the timestamp in the logging
   *
   * @param flag
   */
  void setIncludeTime(const bool flag) {
    m_includeTime = flag;
  }

  /**
   * @brief Set whether to include date
   *
   * @param flag
   */
  void setIncludeDate(const bool flag) {
    m_includeDate = flag;
  }

  /**
   * @brief Get current logging level
   *
   * @return const LoggingLevel&
   */
  const LoggingLevel& currentLevel() const {return m_currentLevel;}

  /**
   * @brief Get global logging level
   *
   * @return const LoggingLevel&
   */
  const LoggingLevel& logLevel() const {return m_logLevel;}

  /**
   * @brief Functor that sets the current logging level
   *
   * @param level Logging level
   * @param file Name of the file where this is called form
   * @param function name of the function
   * @param line Line where this is called from
   * @return Logger&
   */
  Logger& operator() (const LoggingLevel& level,
                      const char* file,
                      int line,
                      const char* function) {
    m_currentLevel = level;
    std::time_t t = std::time(0);
    std::tm* now = std::localtime(&t);

    std::string date("");
    std::string time("");

    const bool addPipe = m_includeDate || m_includeFunction || m_includeLines || m_includeTime;

    if (m_includeTime) {
      time+= " " + formatTime(now->tm_hour) + ":" + formatTime(now->tm_min) + ":" + formatTime(now->tm_sec);
    }

    if (m_includeDate) {
      date += " " + formatTime(now->tm_mday) + "-" + formatTime(now->tm_mon+1) + "-" + std::to_string(now->tm_year+1900);
    }

    if (level <= m_logLevel) {
      m_stream << fancyHeader(level)
               << (m_includeLines ? formatStringFile(file, line, 20) : "")
               << (m_includeFunction ? formatStringFunction(function, 20) : "")
               << (m_includeDate ? date : "")
               << (m_includeTime ? time : "")
               << (addPipe ? " | " : "");
    }
    return *this;
  }

  /**
   * @brief << Operator that allows to use the class an standard streams
   *
   * @tparam T
   * @param l Logger
   * @param message Message to be printed
   * @return Logger&
   */
  template<typename T>
  Logger& operator <<(const T& message) {
    if (this->currentLevel() <= this->logLevel()) {
      std::cout << message;
      return *this;
    } else {
      return *this;
    }
  }
   /**
   * @brief This overload allows std::endl use.
   *
   * @tparam os : Function pointer to the std::endl function.
   * @param
   */
  Logger& operator<< (std::ostream& (*const os)(std::ostream&))
  {
    if (this->currentLevel() <= this->logLevel()) {
      std::cout << os;
      return *this;
    } else {
      return *this;
    }
  }

private:

  /**
   * @brief Construct a new Logger object
   *
   */
  Logger() :
    m_logLevel(LoggingLevel::INFO),
    m_currentLevel(LoggingLevel::INFO),
    m_includeLines(true),
    m_includeFunction(false),
    m_includeDate(true),
    m_includeTime(true) {};

  LoggingLevel m_logLevel;
  LoggingLevel m_currentLevel;
  std::ostream& m_stream = std::cout;
  bool m_includeLines;
  bool m_includeFunction;
  bool m_includeDate;
  bool m_includeTime;

  /**
   * @brief Return nice header based on the current logLevel
   *
   * @param level Log level
   * @return std::string
   */
  static std::string fancyHeader(const LoggingLevel& level) {
    switch (level) {
      case LoggingLevel::ERROR:
        return "\033[0;31m[ ERROR   ]\033[0;0m ";
      case LoggingLevel::WARNING:
        return "\033[0;33m[ WARNING ]\033[0;0m ";
      case LoggingLevel::INFO:
        return "\033[1;32m[ INFO    ]\033[0;0m ";
      case LoggingLevel::DEBUG:
        return "\033[1;36m[ DEBUG   ]\033[0;0m ";
      case LoggingLevel::VERBOSE:
        return "[ VERBOSE ] ";
      default:
        return "";
    }
  }

  /**
   * @brief Format string to always fit in the size
   *
   * @param file input file
   * @param line line number
   * @param max maximum column width
   * @return std::string
   */
  static std::string formatStringFile(const std::string& file, const int line, const std::size_t max) {
    const std::string lineString = std::to_string(line);
    const std::size_t lineSize = lineString.size();
    std::size_t size = (file+":"+lineString).size();
    std::string finalString(file);
    if (size >= max-2) {
      finalString.resize(max-4-lineSize);
      finalString += "...:"+lineString;
    } else {
      finalString += ":" + std::to_string(line) + std::string(max-size, ' ');
    }
    return finalString;
  }

  /**
   * @brief Format string to always fit in the size
   *
   * @param function
   * @param max
   * @return std::string
   */
  static std::string formatStringFunction(std::string function, const std::size_t max) {

    const std::size_t size = function.size();
    if (size >= max-2) {
      function.resize(max - 3);
      function += "...";
    } else {
      function += std::string(max-size, ' ');
    }

    return " " + function;
  }

  /**
   * @brief Turn time to always have 2 digits
   *
   * @param time
   * @return std::string
   */
  static std::string formatTime(int time) {
    if (time > 9) return std::to_string(time);
    return "0" + std::to_string(time);
  }

};

