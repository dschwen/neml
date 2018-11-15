#include "state.h"

namespace neml {

size_t type_size(RealType type) {
  switch(type) {
    case SCALAR:
      return 1;
      break;
    case VECTOR:
      return 3;
      break;
    case MANDEL:
      return 6;
      break;
    case TENSOR:
      return 9;
      break;
  }
}

State::State()
{

}

State::~State()
{
  for (auto it = real_stored_.begin(); it != real_stored_.end(); ++it) {
    if (it->second) delete [] reals_[it->first];
  }
}

std::string State::type()
{
  return "State";
}

std::unique_ptr<NEMLObject> State::initialize(ParameterSet & params)
{
  return make_unique<State>();
}

ParameterSet State::parameters()
{
  ParameterSet pset(State::type());

  return pset;
}

int State::add_real_reference(std::string name, RealType type, 
                              double * location)
{
  real_types_.insert({name, type});
  reals_.insert({name, location});
  real_stored_.insert({name, false});

  return 0;
}

int State::add_real_stored(std::string name, RealType type)
{
  real_types_.insert({name, type}); 
  reals_.insert({name, new double[type_size(type)]});
  real_stored_.insert({name, true});

  return 0;
}

int State::add_state(std::string name)
{
  states_.insert({name, make_unique<State>()});

  return 0;
}


} // namespace neml
