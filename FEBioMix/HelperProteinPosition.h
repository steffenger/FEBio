#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <FECore/FEMaterialPoint.h>
#include "FECore/FESolidDomain.h"
#include "FECore/FEElementSet.h"

double get_distance(FEMaterialPoint& pt_temp, FEMaterialPoint& pt);
double get_distance_shallow(FEMaterialPoint& pt_temp, std::vector<double> point);
std::vector<double> FindNearestMaterialPoint(FEMaterialPoint& pt, FEElementSet* elementSet);
double get_normalized_position(FEMaterialPoint& pt, FEElementSet* elementSetA, FEElementSet* elementSetB);