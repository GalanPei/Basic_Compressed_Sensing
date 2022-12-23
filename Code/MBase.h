#ifndef CODE_MBASE_H
#define CODE_MBASE_H

#include <string>

namespace mmath {

enum class DataType : unsigned {
    INT8 = 0,
    FP16,
    FP32,
};

class MBase {
   public:
    MBase() = default;
    ~MBase() = default;

   private:
    DataType data_type_{DataType::INT8};
};

}  // namespace mmath

#endif