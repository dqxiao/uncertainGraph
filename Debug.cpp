//
//  Debug.cpp
//  uncertainGraph
//
//  Created by dongqingxiao on 10/6/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#include "Debug.hpp"
#include <iostream>
#include <vector>


std::ostream &operator<<(std::ostream &os, std::vector<double> &v)
{
    for(auto &i: v)
        os << v << std::endl;
    
    return os;
}
