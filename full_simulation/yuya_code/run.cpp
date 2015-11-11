//
//  run.cpp
//  AP
//
//  Created by yuya on 2015/09/15.
//  Copyright (c) 2015年 hy. All rights reserved.
//

#include "run.h"


void run() {
    int totalTurn, cellPosition, iE, eE;
    double dr, dtheta, dz;
    double rFieldCheck, zFieldCheck;
    double betaF, betaD;
    clock_t t1, t2;
    t1 = clock();
    /*******************************/
    /*******************************/
    /**** CALCULATION CONDITION ****/
    // fieldMap means which field do you want?
    // Look at getData() in data.cpp
    fieldMap = Kyoto;
    // useToscaField should be always true.
    useToscaField = true;
    
    // essential condition
    // kinetic energy.
    kinetic = 10;// [MeV]
    // initialCondition 0 ~ 5 = {r [m], t [sec], z[m], pr [MeV/c], pth [MeV/c], pz [MeV/c]}
    if (useToscaField) {
        initialCondition[0] = 4.4;
    }
    else {
        initialCondition[0] = sqrt(2 * NatMass(mass) * kinetic + kinetic * kinetic) / ( 300 * 1.0);
    }
    initialCondition[1] = 0; //t
    initialCondition[2] = 0; //z
    initialCondition[3] = 0; //pr
    initialCondition[4] = SIMom(sqrt(2 * NatMass(mass) * kinetic + kinetic * kinetic)); //pth
    initialCondition[5] = 0; //pz
    //sub essential condition
    // dr [m], dz [m] mean deviation from closed orbit.
    dr = 0.001;
    dz = 0;
    // dtheta [deg] -> each step in runge-kutta, converted into rad in calculation.
    dtheta = 0.001;
    // In calculation of tune shift, iE and eE are used.
    iE = kinetic;
    eE = 10;
    // totalTurn -> How many turns do you want?
    totalTurn = 100;
    //sub non-essential condition
    // cellPosition -> later.
    cellPosition = 17;
    // r, zFieldCheck -> youcand check mag field at these points.
    rFieldCheck = 5.1;
    zFieldCheck = 0.02;
    // These two are not necessary now.
    betaF = 0.106 * 2 * 180 / M_PI;
    betaD = 25;
    /****     END     CONDITION ****/
    /*******************************/
    /*******************************/
    
    /************************************************************/
    /************************************************************/
    /************************************************************/
    /*******               MAIN CALCULATION               *******/
    getData();
    getSize(dtheta);
    initialize(initialCondition);
    cout << "start calculating." << endl;
    findClosedOrbit(initialCondition, dtheta, cellPosition);
    drawClosedOrbit(initialCondition, dtheta);
    betatron(initialCondition, dtheta, dr, dz, totalTurn, cellPosition);
    // Checking magnetic field.
    // I'm not sure rocalMagField is correct. 2015/10/19
    rocalMagField(450, 520, 15);
    halfCellMagField(rFieldCheck, zFieldCheck);
    
    // In not beta, the calculation is done by shifting dr, dz from closed orbit separately.
    // In beta, at the same time.
    tuneShiftCalculation(initialCondition, dtheta, dr, dz, iE, eE, totalTurn);
    tuneShiftCalculationBeta(initialCondition, dtheta, dr, dz, iE, eE, totalTurn);
    
    
    /*******             END MAIN CALCULATION             *******/
    /************************************************************/
    /************************************************************/
    /************************************************************/
    
    t2 = clock();
    cout << "finish!!" << endl;
    cout << "calculation time: " << (int)((double)(t2 - t1) / CLOCKS_PER_SEC) / 60 << " m "
    << (int)((double)(t2 - t1) / CLOCKS_PER_SEC) % 60 << "s" << endl;
}