
#include "kernel.h"

//define cubic SPH Kernel
double Wab(double r)
{
    //r /= H;
    std::cout << "r/h = " << r << std::endl;
    if (r < 1.0)
        return (1.0 - (3.0/2.0)*r*r + (3.0/4.0)*r*r*r);//*a2;
    else if (r < 2.0)
        return (0.25*(2.0 - r)*(2.0 - r)*(2.0 - r));//*a2;
    else
        return 0.0;
}

// define gradient of cubic kernel function
double DWab(Point<3,double> & dx, Point<3,double> & DW, double r, bool print)
{
    const double qq=r; ///H;
    double qq2 = qq * qq;
    double fac1 = (c1*qq + d1*qq2)/a2;// /r;
    double b1 = (qq < 1.0)?1.0f:0.0f;
    double wqq = (2.0 - qq);
    double fac2 = c2 * wqq * wqq/a2;// / r;
    double b2 = (qq >= 1.0 && qq < 2.0)?1.0f:0.0f;
    double factor = (b1*fac1 + b2*fac2);
    DW.get(0) = factor * dx.get(0);
    DW.get(1) = factor * dx.get(1);
    DW.get(2) = factor * dx.get(2);
    return factor;
}

// Tensile correction
double Tensile(double r, double rhoa, double rhob, double prs1, double prs2, double temp_wdap)
{
    const double qq=r/H;
    //-Cubic Spline kernel
    double wab;
    if(r>H)
    {
        double wqq1=2.0f-qq;
        double wqq2=wqq1*wqq1;
        wab=a2_4*(wqq2*wqq1);
    }
    else
    {
        double wqq2=qq*qq;
        double wqq3=wqq2*qq;
        wab=a2*(1.0f-1.5f*wqq2+0.75f*wqq3);
    }
    //-Tensile correction.
    double fab=wab*temp_wdap;
    fab*=fab; fab*=fab; //fab=fab^4
    const double tensilp1=(prs1/(rhoa*rhoa))*(prs1>0? 0.01: -0.2);
    const double tensilp2=(prs2/(rhob*rhob))*(prs2>0? 0.01: -0.2);
    return (fab*(tensilp1+tensilp2));
}
