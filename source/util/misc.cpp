/*
misc.cpp
Miscellaneous functions for manipulating data
*/

#include "misc.h"


int find_closest(double x, double array[], int array_length)
{
    /* Traverses the whole array, so O(N) search.  
     * A binary search would be more efficient for long arrays: O(ln N).     
    */
    int index = 0;
    double closest_value = array[0];
    
    for( int i=1; i < array_length; ++i )
    {
        if( std::abs(x - closest_value) > std::abs(x - array[i]) )
        {
            closest_value = array[i];
            index = i;
        }
    }
    
    return index;
}
