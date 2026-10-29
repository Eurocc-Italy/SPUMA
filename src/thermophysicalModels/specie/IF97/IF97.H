/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025-2026 Cineca
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef IF97_H
#define IF97_H

#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace IF97
{

constexpr double R       = 461.526;   // [J/kg/K]
constexpr double Tc      = 647.096;   // [K]
constexpr double pc      = 22.064e6;  // [Pa]
constexpr double rhoc    = 322.0;     // [kg/m3]

#include "../../specie/IF97/reg2meta.H"

// Dynamic viscosity
inline double mu(double rho, double T)
{
    constexpr double muStar = 1.e-6;

    double delta = rho/rhoc;
    double tau = Tc/T;

    double s = 0;
    constexpr double H0[] = {1.67752, 2.20462, 0.6366564, -0.241605};
    for (int i=0; i<4; i++)
        s += H0[i]*pow(tau,i);
    double mu0 = 100/sqrt(tau)/s;


    constexpr double H1[6][7] = {
        { 5.20094E-1, 2.22531E-1, -2.81378E-1, 1.61913E-1, -3.25372E-2, 0.0,         0.0},
        { 8.50895E-2, 9.99115E-1, -9.06851E-1, 2.57399E-1,  0.0,        0.0,         0.0},
        {-1.08374,    1.88797,    -7.72479E-1, 0.0,         0.0,        0.0,         0.0},
        {-2.89555E-1, 1.26613,    -4.89837E-1, 0.0,         6.98452E-2, 0.0,        -4.35673E-3},
        { 0.0,        0.0,        -2.57040E-1, 0.0,         0.0,        8.72102E-3,  0.0},
        { 0.0,        1.20573E-1,  0.0,        0.0,         0.0,        0.0,        -5.93264E-4}
    };

    double s1 = 0;
    for (int i=0; i<6; i++)
    {
        double tau1 = pow(tau - 1, i);
        for (int j=0; j<7; j++)
        {
            if (H1[i][j] != 0.0)
            {
                s1 += H1[i][j]*tau1*pow(delta - 1, j);
            }
        }
    };
    double mu1 = exp(delta*s1);

    return muStar*mu0*mu1;

}

// Interpolation of lambda as in https://srd.nist.gov/JPCRD/jpcrd256.pdf
// Calculation omits the lambda2 correction for region near critical parameters
// Thermal conductivity in [W/mK]

inline double lambda(double rho, double T)
{
    const double tau = T/647.3;
    const double delta = rho/317.7;
    
    constexpr double a[] = {0.0102811, 0.0299621, 0.0156146, -0.00422464};
    double lambda0 = 0;
    for (int i=0; i<4; i++)
        lambda0 += a[i]*pow(tau,i);
    lambda0 *= sqrt(tau);

    constexpr double b0 = -0.397070;
    constexpr double b1 = 0.400302;
    constexpr double b2 = 1.060000;
    constexpr double B1 = -0.171587;
    constexpr double B2 = 2.392190;

    double lambda1 = b0 + b1*delta + b2*exp(B1*(delta + B2)*(delta + B2));
    return lambda0 + lambda1;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace IF97

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

