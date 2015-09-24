/* example.c */
#include "example.h"

/* Compute GCD of two positive integers x and y */
int gcd(int x, int y) {
    int g;
    g = y;
    while (x>0) {
        g = x;
        x = y % x;
        y = g;
    }
    return g;
}