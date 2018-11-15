#ifndef STATE_H
#define STATE_H

#include "objects.h"

#include <map>
#include <string>

namespace neml {

enum RealType { SCALAR, VECTOR, MANDEL, TENSOR };

size_t type_size(RealType type);

class State: public NEMLObject {
 public:
  State();
  virtual ~State();

  static std::string type();
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  static ParameterSet parameters();

  int add_real_reference(std::string name, RealType type, double * location);
  int add_real_stored(std::string name, RealType type);

  int add_state(std::string name);

 private:
  std::map<std::string, double*> reals_;
  std::map<std::string, RealType> real_types_;
  std::map<std::string, bool> real_stored_;
  std::map<std::string, std::unique_ptr<State>> states_;
};

} // namespace neml

#endif // STATE_H
