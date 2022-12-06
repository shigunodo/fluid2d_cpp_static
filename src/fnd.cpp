#include "fnd.h"
#include <algorithm>
#include <cmath>


// definitions of methods for the class Fnd are located



//
//-----------------Reconstruction Scheme------------------
//


// MUSCL: k=1/3 (van Leer, 1977) w/ minmod limitter
// evaluating L/R at cell-boundary between q and qp
void Fnd::MUSCL3(const double& qm, const double& q,
    const double& qp, const double& q2p,
    double& qL, double& qR)
{
    const double k   = 1.0 / 3.0;
    const double b   = (3.0 - k) / (1.0 - k);
    // calculation of L
    double dp  = qp - q;
    double dm  = q - qm;
    double ddp = Fnd::minmod(dp, b * dm);
    double ddm = Fnd::minmod(dm, b * dp);
    qL  = q + 0.25 * ((1.0 - k) * ddm + (1.0 + k) * ddp);
    // calculation of R
    dp  = q2p - qp;
    dm  = qp - q;
    ddp = Fnd::minmod(dp, b * dm);
    ddm = Fnd::minmod(dm, b * dp);
    qR = qp - 0.25 * ((1.0 - k) * ddp + (1.0 + k) * ddm);
}





// MP5: alp=2 (Suresh & Huynh, 1997)
// evaluating L/R at cell-boundary between q and qp
void Fnd::MP5_sub(const double& q2m, const double& qm,
    const double& q, const double& qp,
    const double& q2p, double& qL)
{
    const double alp = 2.0;
    qL = (2.0 * q2m - 13.0 * qm + 47.0 * q + 27.0 * qp - 3.0 * q2p) / 60.0;
    double qMP = q + Fnd::minmod(qp - q, alp * (q - qm));
    if ((qL-q)*(qL-qMP)>1.0e-10)
    {
        double qLL = qL;
        double dm = q2m + q - 2.0 * qm;
        double d = qm + qp - 2.0 * q;
        double dp = q + q2p - 2.0 * qp;
        double dMm = Fnd::minmod4(4.0 * dm - d, 4.0 * d - dm, dm, d);
        double dMp = Fnd::minmod4(4.0 * d - dp, 4.0 * dp - d, d, dp);
        double qUL = q + alp * (q - qm);
        double qMD = 0.5 * (q + qp) - 0.5 * dMp;
        double qLC = q + 0.5 * (q - qm) + 4.0 / 3.0 * dMm;
        double qmin = std::max(std::min({ q,qp,qMD }), std::min({ q,qUL,qLC }));
        double qmax = std::min(std::max({ q,qp,qMD }), std::max({ q,qUL,qLC }));
        qL = Fnd::median(qLL, qmin, qmax);
    }
}