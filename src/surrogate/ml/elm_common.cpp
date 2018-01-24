#include "elm_common.h"

namespace SLFN {

std::ostream& operator<<(std::ostream& os, const ELMActFunc& af) {
    switch (af)
    {
    case ELMActFunc::SIGMOID      : os << "SIGMOID"; break;
    case ELMActFunc::GAUSSIAN     : os << "GAUSSIAN"; break;
    case ELMActFunc::MULTIQUADRIC : os << "MULTIQUADRIC"; break;
    case ELMActFunc::HARDLIMIT    : os << "HARDLIMIT"; break;
    case ELMActFunc::RBF          : os << "RBF"; break;
    case ELMActFunc::TANH         : os << "TANH"; break;
    case ELMActFunc::SIN          : os << "SIN"; break;
    case ELMActFunc::COS          : os << "COS"; break;
    case ELMActFunc::INV_MULTIQUADRIC : os << "INV_MULTIQUADRIC"; break;
    default           : os << "UNKNOWN_ACT_FUNC"; break;
    }
    return os;
}

}
