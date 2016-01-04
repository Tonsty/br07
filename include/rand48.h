#include <stdlib.h>

inline void srand48(long long seed) { srand(seed); }
inline double drand48() { return (double)rand()/RAND_MAX; }