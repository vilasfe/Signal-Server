#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <print>
#include <sys/stat.h>

/*
*  Propagation model test script for signal server
*  Requires gnuplot
*  Compile: gcc -Wall -o test test.cc sui.cc cost.cc ecc33.cc ericsson.cc fspl.cc egli.cc hata.cc -lm
*  Test 850Mhz: ./test 850
*
*  This program is free software; you can redistribute it and/or modify it   *
*  under the terms of the GNU General Public License as published by the     *
*  Free Software Foundation; either version 2 of the License or any later    *
*  version.                                                                  *
*                                                                            *
*  This program is distributed in the hope that it will useful, but WITHOUT  *
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or     *
*  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License     *
*  for more details.                                                         *
*                                                                            *
\****************************************************************************/

#include <string>

#include "cost.hh"
#include "ecc33.hh"
#include "egli.hh"
#include "ericsson.hh"
#include "fspl.hh"
#include "hata.hh"
#include "sui.hh"

extern void point_to_point_ITM(double tht_m, double rht_m, double eps_dielect,
                        double sgm_conductivity, double eno_ns_surfref,
                        double frq_mhz, int radio_climate, int pol,
                        double conf, double rel, double &dbloss, char *strmode,
                        int &errnum);

auto main([[maybe_unused]] int argc, char** argv) -> int
{
double a = 0;
const double f = std::stod(argv[1]);
constexpr double r = 5.0;
constexpr float TxH = 30.0;
constexpr float RxH = 2.0;

mkdir("tests", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

FILE * fh = nullptr;

/*fh = fopen("tests/ITM","w");
for(float d = 0.1; d <= r; d=d+0.2){
    point_to_point_ITM(TxH,RxH,15.0,0.005,301.0,f,5,1,0.5,0.5,a,"",errno);
    fprintf(fh,"%.1f\t%.1f\n",d,a);
}
fclose(fh);
*/
fh = fopen("tests/SUI.1","w");
for(float d = 0.1; d <= r; d=d+0.1){
    a = SUIpathLoss(f, TxH, RxH, d,1);
    std::println(fh,"{:.1f}\t{:.1f}",d,a);
}
fclose(fh);

fh = fopen("tests/COST231.1","w");
for(float d = 0.1; d <= r; d=d+0.1){
    a = COST231pathLoss(f, TxH, RxH, d,1);
    std::println(fh,"{:.1f}\t{:.1f}",d,a);
}
fclose(fh);


fh = fopen("tests/ECC33.1","w");
for(float d = 0.1; d <= r; d=d+0.1){
    a = ECC33pathLoss(f, TxH, RxH, d,1);
    std::println(fh,"{:.1f}\t{:.1f}",d,a);
}
fclose(fh);

fh = fopen("tests/Ericsson9999.1","w");
for(float d = 0.1; d <= r; d=d+0.1){
    a = EricssonpathLoss(f, TxH, RxH, d,1);
    std::println(fh,"{:.1f}\t{:.1f}",d,a);
}
fclose(fh);

fh = fopen("tests/FSPL","w");
for(float d = 0.1; d <= r; d=d+0.1){
    a = FSPLpathLoss(f, d);
    std::println(fh,"{:.1f}\t{:.1f}",d,a);
}
fclose(fh);

fh = fopen("tests/Hata.1","w");
for(float d = 0.1; d <= r; d=d+0.1){
    a = HATApathLoss(f, TxH, RxH, d, 1);
    std::println(fh,"{:.1f}\t{:.1f}",d,a);
}
fclose(fh);

fh = fopen("tests/Egli.VHF-UHF","w");
for(float d = 0.1; d <= r; d=d+0.1){
    a = EgliPathLoss(f, TxH, RxH, d);
    std::println(fh,"{:.1f}\t{:.1f}",d,a);
}
fclose(fh);


fh = fopen("tests/gnuplot.plt","w");
std::println(fh,"set terminal jpeg size 800,600\nset output '{:.0}fMHz_propagation_models.jpg'\nset title '{:.0f}MHz path loss by propagation model - CloudRF.com'",f,f);
std::println(fh,"set key right bottom box\nset style line 1\nset grid\nset xlabel 'Distance KM'\nset ylabel 'Path Loss dB'");
std::print(fh,"plot 'FSPL' with lines, 'Hata.1' with lines, 'COST231.1' with lines,'ECC33.1' with lines,'Ericsson9999.1' with lines,'SUI.1' with lines,'Egli.VHF-UHF' with lines");
fclose(fh);

system("cd tests && gnuplot gnuplot.plt");

return(0);
}

