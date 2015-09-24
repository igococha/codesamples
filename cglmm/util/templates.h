/*
Author: Igor Siveroni

Header file used for the implementation of templated data structures.
Idea taken from stackoverflow "Simulation of templates in C" and Arnold Uthar's blog
with entry "Templates in C.
Adapted to create templates of data structures (METHOD macro).
*/

#ifndef TEMPLATES_H
#define TEMPLATES_H

#define CAT(X,Y) X##_##Y
#define FUNCTION(X,Y) CAT(X,Y)

#define CAT3(X,Y,Z) X##_##Y##_##Z
#define METHOD(X,Y,Z) CAT3(X,Y,Z)


#endif
