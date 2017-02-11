#include <string>
#include <exception>

class ProcessingException: public std::exception {
  private:
    const char * errorMessage;
  public:
    ProcessingException(std::string errorMessage_) {
      errorMessage = errorMessage_.c_str();
      std::cout << errorMessage << std::endl;
    }

    virtual const std::string getMessage() const throw() {
      return errorMessage;
    }
};
