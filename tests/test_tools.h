#ifndef TEST_TOOLS
#define TEST_TOOLS

#define TEST_ASSERT_THROW(condition)                                           \
  {                                                                            \
    if (!(condition)) {                                                        \
      throw std::runtime_error(std::string(__FILE__) + std::string(":") +      \
                               std::to_string(__LINE__) +                      \
                               std::string(" in ") +                           \
                               std::string(__PRETTY_FUNCTION__));              \
    }                                                                          \
  }

#endif // ifndef TEST_TOOLS
