//
//  util.hpp
//  uncertainGraph
//
//  Created by dongqingxiao on 1/18/17.
//  Copyright © 2017 dongqingxiao. All rights reserved.
//

#ifndef util_hpp
#define util_hpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <math.h>
#include <sys/time.h>
#include <armadillo>

using namespace arma;
using namespace std;


int buildSimilarityMatrix(fmat data, fmat& S);

int outputMatrix(fmat matrix, string matrix_filename);

#endif /* util_hpp */
