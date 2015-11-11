//
//  data.h
//  AP
//
//  Created by yuya on 2015/08/02.
//  Copyright (c) 2015年 hy. All rights reserved.
//

#ifndef __AP__data__
#define __AP__data__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include "equation.h"
using namespace std;


void getData();  // put mag field into array from TOSCA file.
void getSize(double);  // some calculation is done in integer. double is converted into integer.
// From current particle position,
// the closest 8 apexes are determined, because of TOSCA's mesh.
//
void getAveFields();

// TOSCA is weird in Coordinate system.
enum coordinateSystem{xyz, rtz, xzy, yxz};
enum map{Kyoto, Kyushu, MERIT, Masaki};

extern map fieldMap;
extern int dataIndex, symmetryNumber;  // dataIndex -> Total number of row of TOSCA data. symmetryNumber -> If KURRI, 12.
// used in getSize()
extern int rDec, thDec, zDec, iRDelta, iThDelta, iZDelta;
extern int dth_decimal;
//
extern int r_size, th_size, z_size;
extern int ir_max, ir_min, iz_max, iz_min;
extern double r_max, r_min, z_max, z_min;
extern bool useToscaField;

// TOSCA data (.LE, .txt, .table) is written in cm and gauss.
// In this code some convertion are done.
extern vector<double> pos[3];   // 0: x,  1: y,  2: z
extern vector<double> mag[3];   // 0: Bx, 1: By, 2: Bz
extern vector<double> dat[6];   // dat -> pos and mag.
extern vector<double> cylin[2]; // 0: r,  1: theta



#endif /* defined(__AP__data__) */
