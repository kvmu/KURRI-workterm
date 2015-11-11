//
//  function.cpp
//  AP
//
//  Created by yuya on 2015/08/28.
//  Copyright (c) 2015年 hy. All rights reserved.
//

#include "function.h"

double rClosedOrbit, rClosedOrbitAt;
double prClosedOrbitAt;
double initialCondition[6];
double closedOrbitLength;
double timeMeasure;


void initialize(double *init) {
    r_var.clear();
    z_var.clear();
    t_var.clear();
    pr.clear();
    pth.clear();
    pz.clear();
    r_var.push_back(init[0]);
    t_var.push_back(init[1]);
    z_var.push_back(init[2]);
    pr.push_back(init[3]);
    pth.push_back(init[4]);
    pz.push_back(init[5]);
    porrection = sqrt(pow(initialCondition[3], 2) + pow(initialCondition[4], 2) + pow(initialCondition[5], 2));
}

void pToXdash(){
    for (int i = 0; i <= dataIndex - 1; i++) {
        pr.at(i) = pr.at(i) / porrection;
        pth.at(i) = pth.at(i) / porrection;
        pz.at(i) = pz.at(i) / porrection;
    }
}

void findClosedOrbit(double *initialCondition, double dth, int rAtTheta) {
    int repeatTime = 0;
    double tmpClosed;
    
    tmpClosed = 0.0;
    initialize(initialCondition);
    if (useToscaField) {
        if (r_var.at(0) > r_max || r_var.at(0) < r_min || z_var.at(0) > z_max || z_var.at(0) < z_min) {
            cout << "r: " << r_var.back() << endl
            << "theta: " << (double)itheta / dth_decimal << endl
            << "z: " << z_var.back() << endl;
            cout << "Divergence" << endl;
            cout << "Magnet Field" << endl;
            cout << "r: " << r_min << " ~ " << r_max << endl;
            cout << "z: " << z_min << " ~ " << z_max << endl;
            exit(0);
        }
    }
repeat:
    for (int ith = 0 * dth_decimal; ith <= 360 / symmetryNumber * dth_decimal; ith += 1) {
        //getAveFieldsで現在の粒子の位置から、実効磁場を求める。
        //dtheta毎に、getAveFieldsをコールする必要がある。
        itheta = ith;
        if (useToscaField) {
            getAveFields();
        }
        
        runge_kutta(drdth, dtdth, dzdth, dprdth, dpthdth, dpzdth, dth * M_PI / 180);
        if (itheta == rAtTheta * dth_decimal && rAtTheta != 0) {
            rClosedOrbitAt = r_var.back();
            prClosedOrbitAt = pr.back();
        }
        
        //Divergence action
        if (useToscaField) {
            if (r_var.back() > r_max || r_var.back() < r_min || z_var.back() > z_max || z_var.back() < z_min) {
                cout << "r: " << r_var.back() << endl
                << "theta: " << (double)itheta / dth_decimal << endl
                << "z: " << z_var.back() << endl;
                cout << "Divergence" << endl;
                cout << "Magnet Field" << endl;
                cout << "r: " << r_min << " ~ " << r_max << endl;
                cout << "z: " << z_min << " ~ " << z_max << endl;
                exit(0);
            }
        }
    }
    
    if (useToscaField) {
        if (fabs(r_var.at(0) - r_var.back()) > 1e-10 || NatMom(fabs(pth.at(0) - pth.back())) > 1e-10) {
            if (repeatTime > 100) {
                cout << "Error may large." << endl;
                cout << "pr at    0 deg: " << setprecision(10) << NatMom(pth.at(0)) << endl;
                cout << "pr at cell deg: " << setprecision(10) << NatMom(pth.back()) << endl;
                exit(0);
            }
            tmpClosed = (r_var.at(0) + r_var.back()) / 2;
            initialize(initialCondition);
            r_var.at(0) = tmpClosed;
            repeatTime++;
            goto repeat;
        }
        else {
            tmpClosed = r_var.back();
            cout << "closed orbit found." << endl
            << "initial r: " << fixed << setprecision(10) << tmpClosed << endl;
        
            rClosedOrbit = tmpClosed;
        }
    }
    else if (!useToscaField) {
        // Weak focusing
        rClosedOrbit = 0.456557;
    }

}

void getClosedOrbitLength(double *initialCondition, double dth) {
    double length = 0;
    double dtheta = M_PI / (180 * dth_decimal);
    closedOrbitLength = 0;
    initialCondition[4] = SIMom(sqrt(2 * NatMass(mass) * kinetic + kinetic * kinetic)); //pth
    initialize(initialCondition);
    r_var.at(0) = rClosedOrbit;
    for (int ith = 0 * dth_decimal; ith <= 360 * dth_decimal; ith += 1) {
        itheta = ith;
        if (useToscaField) {
            getAveFields();
        }
        runge_kutta(drdth, dtdth, dzdth, dprdth, dpthdth, dpzdth, dth * M_PI / 180);
        timeMeasure += t_var.back();
        length += r_var.back() * dtheta;
    }
    
    closedOrbitLength = length;
    
    
}


void drawClosedOrbit(double *initialCondition, double dtheta) {
    ofstream closedOrbit("closedOrbit.csv");
    ofstream excursion("excursion.csv");
    initialize(initialCondition);
    r_var.at(0) = rClosedOrbit;
    z_var.at(0) = 0;
    closedOrbit << "x [cm]" << "," << "y [cm]" << endl;
    for (int ith = 0; ith <= 360 * dth_decimal; ith += 1) {
        itheta = ith;
        if (useToscaField) {
            getAveFields();
        }
        if (ith % 10 == 0) {
            closedOrbit << r_var.back() * cos((double)ith / dth_decimal * M_PI / 180) * 100 << ","
            << r_var.back() * sin((double)ith / dth_decimal * M_PI / 180) * 100 << endl;
            
            if (ith <= 360 / symmetryNumber * dth_decimal) {
                excursion << (double)ith / dth_decimal << "," << r_var.back() * 100 << endl;
            }
        }
        
        runge_kutta(drdth, dtdth, dzdth, dprdth, dpthdth, dpzdth, dtheta * M_PI / 180);
    }
    cout << "excursion " << 100 * (*max_element(r_var.begin(), r_var.end()) - *min_element(r_var.begin(), r_var.end())) << endl;
}
void BL(double betaF, double betaD, double dtheta) {
    ofstream BLF("BLF.csv");
    ofstream BLD("BLD.csv");
    int bF, bD, bS, bCell;
    int r0 = (int)(r_min * 10 + 1e-7);
    int r1 = (int)(r_max * 10 + 1e-7);
    double sum;
    
    bCell = 360 / symmetryNumber * dth_decimal;
    bF = (int)(betaF * dth_decimal + 1e-7);
    bD = (int)(betaD * dth_decimal + 1e-7);
    bS = (bCell - bD - bF - bD) / 2;
    sum = 0;
    
    for (int r = r0; r <= r1; r++) {
        for (int ith = bS; ith <= bS + bD; ith++) {
            itheta = ith;
            r_var.push_back((double)r / 10);
            z_var.push_back(0);
            getAveFields();
            runge_kutta(drdth, dtdth, dzdth, dprdth, dpthdth, dpzdth, dtheta * M_PI / 180);
            sum += fluz * r_var.back() * dtheta * M_PI / 180;
        }
        BLF << r << "," << sum << endl;
    }
    
    for (int r = r0; r <= r1; r++) {
        for (int ith = bS + bD; ith <= bS + bD + bF; ith++) {
            itheta = ith;
            r_var.push_back((double)r / 10);
            z_var.push_back(0);
            getAveFields();
            runge_kutta(drdth, dtdth, dzdth, dprdth, dpthdth, dpzdth, dtheta * M_PI / 180);
            sum += fluz * r_var.back() * dtheta * M_PI / 180;
        }
        BLD << r << "," << sum << endl;
    }
    
}

void rocalMagField(double r0, double r1, double th) {
    ofstream field("rocalField.csv");
    /*
    for (int r = r0; r <= r1; r++) {
        itheta = th * dth_decimal;
        r_var.push_back((double)r / 100);
        z_var.push_back(0);
        getAveFields();
        field << r_var.back() << "," << z_var.back() << "," << (double)itheta / dth_decimal << ","
        << flux << "," << fluy << "," << fluz << endl;
    }*/
    for (int i = r_size * 20 + 200; i <= dataIndex - 1; i += r_size) {
        field << cylin[0][i] << "," << cylin[1][i] << "," << pos[2][i] << ","
        << mag[0][i] << "," << mag[1][i] << "," << mag[2][i] << endl;
    }
}
void halfCellMagField(double r, double z) {
    ofstream field ("halfCellField.csv");
    
    for (int ith = 0; ith <= 180 / symmetryNumber * dth_decimal; ith += 1) {
        if (ith % 10 != 0) continue;
        itheta = ith;
        r_var.push_back(r);
        z_var.push_back(z);
        getAveFields();
        field << r_var.back() << "," << z_var.back() << "," << (double)itheta / dth_decimal << ","
        << flux << "," << fluy << "," << fluz << endl;
    }
    r_var.clear();
    z_var.clear();
}

void betatron(double * initial_condition, double dtheta, double dr, double dz, int totalTurn, int cellPosition) {
    int tmp_theta, thPeriod;
    ofstream lattice("lattice.csv");
    initialize(initial_condition);
    r_var.at(0) = rClosedOrbit + dr;
    z_var.at(0) = dz;
    thPeriod = cellPosition * dth_decimal;
    
    lattice << "theta [deg]" << "," << "r [m]" << "," << "z [m]" << ","
    << "pr [MeV/c]" << "," << "pth [Mev/c]" << "," << "pz [MeV/c]" << "," << "pr / pth" << endl;
    for (int ith = 0; ith <= 360 * totalTurn * dth_decimal; ith += 1) {
        //getAveFieldsで現在の粒子の位置から、実効磁場を求める。
        //dtheta毎に、getAveFieldsをコールする必要がある。
        itheta = ith;
        tmp_theta = ith;

        if (useToscaField) {
            getAveFields();
        }
        if (ith % (360 * dth_decimal) == 0 && ith != 0) {
            cout << ith / (360 * dth_decimal) << endl;
        }
        
        if (tmp_theta == thPeriod && (tmp_theta != 0)) {
            lattice << tmp_theta / dth_decimal << ","
            << fixed << setprecision(20)
            << r_var.back() - rClosedOrbitAt << ","
            << z_var.back() << ","
            << NatMom(pr.back()) - NatMom(prClosedOrbitAt) << ","
            << NatMom(pth.back()) << ","
            << NatMom(pz.back()) << ","
            << NatMom(pr.back()) / NatMom(pth.back()) << endl;
            
            thPeriod += 360 / symmetryNumber * dth_decimal;
        }
        
        runge_kutta(drdth, dtdth, dzdth, dprdth, dpthdth, dpzdth, dtheta * M_PI / 180);
    }
}

void tuneShiftCalculation(double *initialCondition, double dtheta, double dr, double dz, int iE, int eE, int totalTurn) {
    double tuneHorizontal, tuneVertical;
    ofstream tuneShift("tuneShift.csv");
    
    tuneShift << "Qh" << "," << "Qv" << endl;
    for (int ikinetic = iE; ikinetic <= eE; ikinetic++) {
        kinetic = ikinetic;
        initialCondition[4] = SIMom(sqrt(2 * NatMass(mass) * kinetic + kinetic * kinetic)); //pth
        porrection = sqrt(pow(initialCondition[3], 2) + pow(initialCondition[4], 2) + pow(initialCondition[5], 2));
        findClosedOrbit(initialCondition, dtheta, 0);
        tuneHorizontal = Hune(initialCondition, dtheta, dr, totalTurn);
        tuneVertical = Vune(initialCondition, dtheta, dz, totalTurn);
        tuneShift << tuneHorizontal << "," << tuneVertical << endl;
    }
}


void tuneShiftCalculationBeta(double *initialCondition, double dtheta, double dr, double dz, int iE, int eE, int totalTurn) {
    double tuneHorizontal, tuneVertical;
    ofstream tuneShift("tuneShiftBeta.csv");
    
    tuneShift << "Qh" << "," << "Qv" << endl;
    for (int ikinetic = iE; ikinetic <= eE; ikinetic++) {
        kinetic = ikinetic;
        initialCondition[4] = SIMom(sqrt(2 * NatMass(mass) * kinetic + kinetic * kinetic)); //pth
        porrection = sqrt(pow(initialCondition[3], 2) + pow(initialCondition[4], 2) + pow(initialCondition[5], 2));
        findClosedOrbit(initialCondition, dtheta, 0);
        Tune(initialCondition, dtheta, dr, dz, totalTurn, tuneHorizontal, tuneVertical);
        tuneShift << tuneHorizontal << "," << tuneVertical << endl;
    }
}

void Tune(double * initial_condition, double dtheta, double dr, double dz, int totalTurn, double &htune, double &vtune) {
    int tmp_theta, advCell, thPeriod;
    double summuh, summuv;
    initialize(initial_condition);
    r_var.at(0) = rClosedOrbit + dr;
    z_var.at(0) = dz;
    summuh = 0;
    summuv = 0;
    htune = 0;
    vtune = 0;
    thPeriod = 360 / symmetryNumber * dth_decimal;
    cout << "calculating Tune at " << (int)(kinetic + 1e-4) << " MeV ...." << endl;
    for (int ith = 0; ith <= 360 * totalTurn * dth_decimal; ith += 1) {
        //getAveFieldsで現在の粒子の位置から、実効磁場を求める。
        //dtheta毎に、getAveFieldsをコールする必要がある。
        itheta = ith;
        tmp_theta = ith;
        
        if (useToscaField) {
            getAveFields();
        }
        
        if (ith % (360 * dth_decimal) == 0 && ith != 0) {
            cout << ith / (360 * dth_decimal) << endl;
        }
        
        if ((tmp_theta % (360 / symmetryNumber * dth_decimal)) == 0 && (tmp_theta != 0)) {
            summuh += acos(
                          ((rClosedOrbit - r_var.at(ith)) * (rClosedOrbit - r_var.at(ith - thPeriod)) + pr.at(ith) * pr.at(ith - thPeriod))
                          / (sqrt((rClosedOrbit - r_var.at(ith)) * (rClosedOrbit - r_var.at(ith)) + pr.at(ith) * pr.at(ith))
                             * sqrt((rClosedOrbit - r_var.at(ith - thPeriod)) * (rClosedOrbit - r_var.at(ith - thPeriod)) + pr.at(ith - thPeriod) * pr.at(ith - thPeriod)))
                          );
            summuv += acos(
                          ((z_var.at(ith)) * (z_var.at(ith - thPeriod)) + pz.at(ith) * pz.at(ith - thPeriod))
                          / (sqrt(z_var.at(ith) * z_var.at(ith) + pz.at(ith) * pz.at(ith))
                             * sqrt(z_var.at(ith - thPeriod) * z_var.at(ith - thPeriod) + pz.at(ith - thPeriod) * pz.at(ith - thPeriod)))
                          );
        }
        
        runge_kutta(drdth, dtdth, dzdth, dprdth, dpthdth, dpzdth, dtheta * M_PI / 180);
    }
    advCell = symmetryNumber * totalTurn;
    htune = summuh * symmetryNumber / (advCell * 2 * M_PI);
    vtune = summuv * symmetryNumber / (advCell * 2 * M_PI);
}

double Hune(double * initial_condition, double dtheta, double dr, int totalTurn) {
    int tmp_theta, advCell, thPeriod;
    double summu, HTune;
    initialize(initial_condition);
    r_var.at(0) = rClosedOrbit + dr;
    summu = 0;
    thPeriod = 360 / symmetryNumber * dth_decimal;
    cout << "calculating Hune at " << (int)(kinetic + 1e-4) << " MeV ...." << endl;
    for (int ith = 0; ith <= 360 * totalTurn * dth_decimal; ith += 1) {
        //getAveFieldsで現在の粒子の位置から、実効磁場を求める。
        //dtheta毎に、getAveFieldsをコールする必要がある。
        itheta = ith;
        tmp_theta = ith;
        
        if (useToscaField) {
            getAveFields();
        }
        if (ith % (360 * dth_decimal) == 0 && ith != 0) {
            cout << ith / (360 * dth_decimal) << endl;
        }
        
        if ((tmp_theta % (360 / symmetryNumber * dth_decimal)) == 0 && (tmp_theta != 0)) {
            summu += acos(
                          ((rClosedOrbit - r_var.at(ith)) * (rClosedOrbit - r_var.at(ith - thPeriod)) + pr.at(ith) * pr.at(ith - thPeriod))
                          / (sqrt((rClosedOrbit - r_var.at(ith)) * (rClosedOrbit - r_var.at(ith)) + pr.at(ith) * pr.at(ith))
                             * sqrt((rClosedOrbit - r_var.at(ith - thPeriod)) * (rClosedOrbit - r_var.at(ith - thPeriod)) + pr.at(ith - thPeriod) * pr.at(ith - thPeriod)))
                          );
        }
        
        runge_kutta(drdth, dtdth, dzdth, dprdth, dpthdth, dpzdth, dtheta * M_PI / 180);
    }
    advCell = symmetryNumber * totalTurn;
    HTune = summu * symmetryNumber / (advCell * 2 * M_PI);
    
    return HTune;
}

double Vune(double * initial_condition, double dtheta, double dz, int totalTurn) {
    int tmp_theta, advCell, thPeriod;
    double summu, VTune;
    initialize(initial_condition);
    r_var.at(0) = rClosedOrbit;
    z_var.at(0) = dz;
    summu = 0;
    thPeriod = 360 / symmetryNumber * dth_decimal;
    cout << "calculating Vune at " << (int)(kinetic + 1e-4) << " MeV ...." << endl;
    for (int ith = 0; ith <= 360 * totalTurn * dth_decimal; ith += 1) {
        //getAveFieldsで現在の粒子の位置から、実効磁場を求める。
        //dtheta毎に、getAveFieldsをコールする必要がある。
        itheta = ith;
        tmp_theta = ith;
        
        if (useToscaField) {
            getAveFields();
        }
        if (ith % (360 * dth_decimal) == 0 && ith != 0) {
            cout << ith / (360 * dth_decimal) << endl;
        }
        
        if ((tmp_theta % (360 / symmetryNumber * dth_decimal)) == 0 && (tmp_theta != 0)) {
            summu += acos(
                          ((z_var.at(ith)) * (z_var.at(ith - thPeriod)) + pz.at(ith) * pz.at(ith - thPeriod))
                          / (sqrt(z_var.at(ith) * z_var.at(ith) + pz.at(ith) * pz.at(ith))
                             * sqrt(z_var.at(ith - thPeriod) * z_var.at(ith - thPeriod) + pz.at(ith - thPeriod) * pz.at(ith - thPeriod)))
                          );
        }
        
        runge_kutta(drdth, dtdth, dzdth, dprdth, dpthdth, dpzdth, dtheta * M_PI / 180);
    }
    advCell = symmetryNumber * totalTurn;
    VTune = summu * symmetryNumber / (advCell * 2 * M_PI);
    
    return VTune;
}