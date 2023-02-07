#include "atomic_number.h"
#include <math.h>

/*
 * Get the atomic number by selecting the atom whose atomic mass is the closest
 * to the atomic mass that is used. This has the advantage of working for all
 * isotopes with the possible exception of tritium (which is closer to He).
 */
int get_atomic_number(REAL_T mass) {
   int i, atomic_num = 0;
   REAL_T min_diff = fabs(mass);
   REAL_T diff;

   for (i = 0; i < NUM_ATOMS; i++) {
      diff = fabs(mass - ATOMIC_MASS[i]);
      if (diff < min_diff) {
         min_diff = diff;
         atomic_num = i;
      }
   }

   return atomic_num;
}
