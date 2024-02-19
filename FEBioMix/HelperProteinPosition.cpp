#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include "HelperProteinPosition.h"

double get_distance(FEMaterialPoint& pt_temp, FEMaterialPoint& pt){
return std::sqrt(std::pow(pt_temp.m_r0.x - pt.m_r0.x, 2) + std::pow(pt_temp.m_r0.y - pt.m_r0.y, 2));
}

double get_distance_shallow(FEMaterialPoint& pt_temp, std::vector<double> point){
return std::sqrt(std::pow(pt_temp.m_r0.x - point[0], 2) + std::pow(pt_temp.m_r0.y - point[1], 2));
}

std::vector<double> FindNearestMaterialPoint(FEMaterialPoint& pt, FEElementSet* elementSet){
std::vector<double> nearestMaterialPoint_shallowCopy = {0, 0};
double minDistance = std::numeric_limits<double>::max();
    
for (int i = 0; i < elementSet->Elements(); i++) {
        FEElement &element = elementSet->Element(i);
        for (int j = 0; j < element.GaussPoints(); j++) {

            FEMaterialPoint &materialPoint = *element.GetMaterialPoint(j);
            double distance = get_distance(pt, materialPoint);
            //std::cout << distance << std::endl;
            if (distance < minDistance) {
               minDistance = distance;
               //std::cout << "mD " << minDistance << std::endl;
               nearestMaterialPoint_shallowCopy[0] = materialPoint.m_r0.x;
               nearestMaterialPoint_shallowCopy[1] = materialPoint.m_r0.y;
            }
            
        }
       }

    return nearestMaterialPoint_shallowCopy;
}


double get_normalized_position(FEMaterialPoint& pt, FEElementSet* elementSetCV, FEElementSet* elementSetPF){
   double protein;
   std::vector<double> CV = FindNearestMaterialPoint(pt, elementSetCV);
   std::vector<double> PF = FindNearestMaterialPoint(pt, elementSetPF);
   double distance1 = get_distance_shallow(pt, CV);
   double distance2 = get_distance_shallow(pt, PF);
   protein = distance2/(distance2+distance1);
   return protein;
}