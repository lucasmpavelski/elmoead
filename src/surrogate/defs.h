#pragma once

#include <cmath>
using std::abs;
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <type_traits>
using std::size_t;
#include <limits>
using std::numeric_limits;
#include <algorithm>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <memory>
#include <iterator>
#include <iterator>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>

#include "Eigen/Core"

/* common used aliases */
template<typename T,int no_rows=Eigen::Dynamic,int no_cols=Eigen::Dynamic,int mj=Eigen::RowMajor>
using ematrix = Eigen::Matrix<T,no_rows,no_cols,mj>;
using ematrixd = ematrix<double>;
using ematrixi = ematrix<int>;
using ematrixb = ematrix<bool>;

template<typename T,int no_el=Eigen::Dynamic>
using evector = Eigen::Matrix<T,no_el,1,Eigen::ColMajor>;
using evectord = evector<double>;
using evectori = evector<int>;
using evectorb = evector<bool>;

template<typename T>
using svector = std::vector<T>;
using svectord = svector<double>;
using svectori = svector<int>;
using svectorb = svector<bool>;
