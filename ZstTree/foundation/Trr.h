#ifndef _FOUNDATION_TRR_H
#define _FOUNDATION_TRR_H

#include "BiStates.h"
#include "BaseDefine.h"
#include "MyInterval.h"
#include "Coordinate.h"

template <typename NN, bool DOCHECK = true>
class F_Trr//todo 为什么感觉这里只有正的矩形，没有斜的
{
public:
    void MakeDiamond ( const F_Point<NN, DOCHECK> p, NN r ) {//todo 为什么这么计算
        NN tval;

        tval = p.x - p.y;
        xlow = tval - r;
        xhi  = tval + r;

        tval = p.x + p.y;
        ylow = tval - r;
        yhi = tval + r;
        
        z = p.z;
    }
    
    F_Trr () {
        xlow = ylow = 1;
        xhi  = yhi  = 0;
        z = 0;
    }
    F_Trr ( NN x1, NN x2, NN y1, NN y2, int z_) {
        xlow = x1;
        xhi = x2;
        ylow = y1;
        yhi = y2;
        z = z_;
    }
    F_Trr ( const F_Point<NN, DOCHECK> p, NN r ) {
        MakeDiamond ( p, r ) ;
    }

    void Enclose ( const F_Trr<NN, DOCHECK> b ) {
        if ( IsEmpty() ) {
            xlow = b.xlow ;
            ylow = b.ylow ;
            xhi  = b.xhi  ;
            yhi  = b.yhi  ;
        } else {
            xlow = tMIN (xlow, b.xlow) ;
            ylow = tMIN (ylow, b.ylow) ;
            xhi  = tMAX (xhi , b.xhi ) ;
            yhi  = tMAX (yhi , b.yhi ) ;
        }
    }
    bool IsEmpty ( ) const {
        F_Interval<NN, DOCHECK> x_interval (xlow, xhi);
        F_Interval<NN, DOCHECK> y_interval (ylow, yhi);
        return x_interval.IsEmpty() || y_interval.IsEmpty() ;
    }
    NN Width ( TwoStates i ) const {
        if ( i == X ) { 
            F_Interval<NN, DOCHECK> x_interval (xlow, xhi);
            if ( x_interval.IsEmpty() ) {
                return 0;
            } else {
                return xhi - xlow;
            }
        } else if ( i == Y ) {
            F_Interval<NN, DOCHECK> y_interval (ylow, yhi);
            if ( y_interval.IsEmpty() ) {
                return 0;
            } else {
                return yhi - ylow;
            }
        } else { exit(0); }
    }
    
    void ErrCorrection ( ) {//todo 这种处理是为了避免FUZZ吗，即将差小于FUZZ的情况直接改成等于？
        NN a = xlow ;
        NN b = xhi  ;
        if ( a == b ) {
            xlow = xhi = (a+b)/2 ;
        }
        a = ylow ;
        b = yhi  ;
        if ( a == b ) {
            ylow = yhi = (a+b)/2 ;
        }
    }
    void SelfDOCHECK ( ) {
        ErrCorrection () ;
        assert ( !IsEmpty() ) ;
    }

    void MakeIntersect ( F_Trr<NN, DOCHECK> trr1, F_Trr<NN, DOCHECK> trr2 ) {
        xlow = tMAX ( trr1.xlow, trr2.xlow ) ;
        xhi  = tMIN ( trr1.xhi , trr2.xhi  ) ;
        ylow = tMAX ( trr1.ylow, trr2.ylow ) ;
        yhi  = tMIN ( trr1.yhi , trr2.yhi  ) ;
        if (trr1.z == trr2.z) {
            z = trr1.z;
        } else {
            z = 0;
        }
        if (DOCHECK) 
            SelfDOCHECK () ;
    }

public:
    NN xlow, xhi, ylow, yhi;
    int z;
};

#endif
