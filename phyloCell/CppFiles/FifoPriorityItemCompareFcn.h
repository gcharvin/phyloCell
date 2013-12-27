// Copyright 2006 The MathWorks, Inc.

#ifndef FifoPriorityItemCompareFcn_h
#define FifoPriorityItemCompareFcn_h

#include "FifoPriorityItem.h"
#include "FifoPriorityQueue.h"

#include "mex.h"

/** FifoPriorityItemCompareFcn

 *  A function object used to control the sort order in the class
 *  FifoPriorityQueue.  The () overload first sorts based on priority,
 *  either highest priority first or lowest priority first.  If
 *  priorites are equal, it sorts based on order, with older items
 *  sorted first.
 *
 *  @see FifoPriorityQueue, FifoPriorityItem
 */
template <typename DATA_T, typename PRIORITY_T>
class FifoPriorityItemCompareFcn
{
  public:
    enum CompareMode {HighestPriorityFirst, LowestPriorityFirst};

    /** FifoPriorityItemCompareFcn
     *  Constructor for comparison function. Default construction uses
     *  HighestPriorityFirst.
     *
     *  @param m    CompareMode, either HighestPriorityFirst or LowestPriorityFirst
     */
    FifoPriorityItemCompareFcn (CompareMode m=HighestPriorityFirst) : fMode(m) {
    }

    /** operator()
     *  Overload implementing "less-than".
     *
     *  @param p1   const FifoPriorityItem<DATA_T, PRIORITY_T>&
     *  @param p2   const FifoPriorityItem<DATA_T, PRIORITY_T>&
     *
     *  @return     bool, true if p1 is "less than" p2.
     */
    bool operator() (const FifoPriorityItem<DATA_T, PRIORITY_T>& p1,
                     const FifoPriorityItem<DATA_T, PRIORITY_T>& p2) const {

        if (getCompareMode() == HighestPriorityFirst) {
            if (p1.getPriority() < p2.getPriority()) {
                return true;
            }
            else if (p1.getPriority() > p2.getPriority()) {
                return false;
            }
        } 
        else {  
            // LowestPriorityFirst
            if (p1.getPriority() > p2.getPriority()) {
                return true;
            }
            else if (p1.getPriority() < p2.getPriority()) {
                return false;
            }
        } 

        // If we get here, the priorities are equal.  Use order
        // as a tie-breaker.

        return (p1.getOrder() > p2.getOrder());
    }

  private:
    CompareMode fMode;   /**< Comparison mode */

    /** getCompareMode
     *  Get the compare mode.
     *
     *  @return   CompareMode
     */
    CompareMode getCompareMode(void) const {
        return fMode;
    }

};

#endif
