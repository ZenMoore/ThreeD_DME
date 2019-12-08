#ifndef _FOUNDATION_MYINTERVAL_H
#define _FOUNDATION_MYINTERVAL_H

#include "BiStates.h"
#include "BaseDefine.h"

template <class NN, bool DOCHECK = true>
class F_Interval {

private:
    NN m_bound[2];

public: 

    /* constructor */
    F_Interval () { m_bound[LOW] = 1; m_bound[HIGH] = 0; }

    explicit F_Interval (NN a) { m_bound[LOW] = a; m_bound[HIGH] = a; }

    F_Interval (NN a, NN b) {
        m_bound[LOW] = a;
        m_bound[HIGH] = b;
    }

    F_Interval (const F_Interval<NN, DOCHECK> &x) {
        m_bound[LOW] = x.m_bound[LOW];
        m_bound[HIGH] = x.m_bound[HIGH];
    }

    /* member function */
    bool IsEmpty () const {
        return m_bound[LOW] > m_bound[HIGH];
    }
    
    bool IsPoint () const {
        return m_bound[LOW] == m_bound[HIGH];
    }

    void Enclose (NN n) {
        if ( IsEmpty() ) {
            m_bound[LOW]  = n ;
            m_bound[HIGH] = n ;
        } else {
            m_bound[LOW]  = tMIN( m_bound[LOW], n );
            m_bound[HIGH] = tMAX( m_bound[HIGH], n );
        }
    }

    void Enclose ( const F_Interval<NN, DOCHECK>& p ) {
        if ( !p.IsEmpty() ) {
            Enclose ( p.m_bound[LOW] ) ;
            Enclose ( p.m_bound[HIGH] ) ;
        }
    }
    
    bool IsEnclose (NN n) const {
        return n >= m_bound[LOW] && n <= m_bound[HIGH];
    }
    
    bool IsEnclose (const F_Interval<NN, DOCHECK>& p) const {
        return p.m_bound[LOW] >= m_bound[LOW] && p.m_bound[HIGH] <= m_bound[HIGH];
    }
    
    NN Width () const {
        if (IsEmpty()) return 0;
        return m_bound[HIGH] - m_bound[LOW];
    }
};

#endif
