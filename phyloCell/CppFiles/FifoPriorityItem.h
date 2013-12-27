// Copyright 2006 The MathWorks, Inc.

#ifndef FifoPriorityItem_h
#define FifoPriorityItem_h

#include "mex.h"

/** FifoPriorityItem
 *  Groups a data value, a priority value, and an order value.
 *
 *  @see FifoPriorityQueue, FifoPriorityItemCompareFcn
 */
template <typename DATA_T, typename PRIORITY_T>
class FifoPriorityItem
{
  public:
   /** FifoPriorityItem
     *  Class constructor.
     *
     *  @param data        DATA_T
     *  @param priority    PRIORITY_T
     *  @param order       uint64_T
     *
     *  @see FifoPriorityQueue, FifoPriorityItemCompareFcn
     */
    FifoPriorityItem(DATA_T data, PRIORITY_T priority, uint64_T order) {
        fData     = data;
        fPriority = priority;
        fOrder    = order;
    }

    /** getData
     *  Returns the item's data value.
     *
     *  @return   DATA_T
     */
    DATA_T getData(void) const {
        return fData;
    }

    /** getPriority
     *  Returns the item's priority value.
     *
     *  @return   PRIORITY_T
     */
    PRIORITY_T getPriority(void) const {
        return fPriority;
    }

    /** getOrder
     *  Returns the item's order.
     *
     *  @return   uint64_T
     */
    uint64_T getOrder(void) const {
        return fOrder;
    }

  private:
    DATA_T     fData;       /**< item's data     */
    PRIORITY_T fPriority;   /**< item's priority */
    uint64_T   fOrder;      /**< item's order    */
};

#endif
