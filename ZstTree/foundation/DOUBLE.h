#ifndef _FOUNDATION_DOUBLE_H
#define _FOUNDATION_DOUBLE_H

#include <cmath>
#include "BaseDefine.h"

class DOUBLE
{
public:    
    double value;

public:
    DOUBLE() {
        value = 0;
    }
    DOUBLE(double v) {
        value = v;
    }

    DOUBLE operator+(DOUBLE b) const {
        return DOUBLE(value + b.value);
    }
    void operator+=(double b) {
        value += b;
    }
    void operator+=(DOUBLE b) {
        value += b.value;
    }
    DOUBLE operator-(DOUBLE b) const {
        return DOUBLE(value - b.value);
    }
    DOUBLE operator*(DOUBLE b) const {
        return DOUBLE(value * b.value);
    }
    void operator*=(DOUBLE b) {
        value *= b.value;
    }
    void operator*=(double b) {
        value *= b;
    }
    DOUBLE operator*(double b) {
        return DOUBLE(value * b);
    }
    DOUBLE operator/(DOUBLE b) const {
        return DOUBLE(value / b.value);
    }
    bool operator<(DOUBLE b) const {
        return value < b.value;
    }
    bool operator<=(DOUBLE b) const {
        return value <= b.value;
    }
    bool operator>=(double b) const {
        return value >= b;
    }
    bool operator>(DOUBLE b) const {
        return value > b.value;
    }
    bool operator>=(DOUBLE b) const {
        return value >= b.value;
    }
    bool operator==(DOUBLE b) const {
        return ABS(value - b.value) <= FUZZ;
    }
    bool operator!=(DOUBLE b) const {
        return (ABS(value - b.value) > FUZZ);
    }
    DOUBLE _ABS() {
        return value >= 0 ? DOUBLE(value) : DOUBLE(-value);
    }
    DOUBLE _SQRT() {
        return DOUBLE(sqrt(value));
    }
    int _ROUND() {
        return floor(value + 0.5);
    }
};

#endif