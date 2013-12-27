// FifoPriorityQueue.h
// Copyright 2006 The MathWorks, Inc.

#ifndef FifoPriorityQueue_H
#define FifoPriorityQueue_H

#include "FifoPriorityItem.h"
#include "FifoPriorityItemCompareFcn.h"
#include <queue>

#include "mex.h"

/** FifoPriorityQueue
 *  Like the standard library's priority queue, but with the following
 *  additional behavior:
 *  Clients push a pairs of values onto the queue: a data value and a priority value.
 *  The priority value determines the order in which items are retrieved from the queue.
 *  Clients specify at construction time whether higher-priority or lower-priority items
 *  are returned first.
 *
 *  Unlike the standard C++ library's priority_queue, this class has FIFO behavior.
 *  If two items have the same priority, the one pushed onto the queue first will be
 *  returned first.  In the priority_queue class, the retrieval order is undefined
 *  for equal priorities.
 *
 *  @see FifoPriorityItem, FifoPriorityItemCompareFcn
 */
template <typename DATA_T, typename PRIORITY_T>
class FifoPriorityQueue
{
  public:
    /** FifoPriorityQueue
     *  Class constructor.
	 *
	 * @param mode enum FifoPriorityItemCompareFcn<DATA_T, PRIORITY_T>::CompareMode
	 *             HighestPriorityFirst specifies that when you get item info from the
	 *             queue using topData() or topPriority(), items with the highest
	 *             priority are returned first.
	 *             LowestPriorityFirst specifies that items with the lowest priority
	 *             are returned first.
     */
    FifoPriorityQueue(enum FifoPriorityItemCompareFcn<DATA_T, PRIORITY_T>::CompareMode mode) {
        FifoPriorityItemCompareFcn<DATA_T, PRIORITY_T> compareFcn(mode);
        queue_T queue(compareFcn);
        setPriorityQueue(queue);
        setCount(0);
    }

    /** push
     *  Push a data value and an associated priority value onto the queue.
     *
     *  @param data         DATA_T
     *  @param priority     PRIORITY_T
     */
    void push(DATA_T data, PRIORITY_T priority) {
        getPriorityQueue().push(FifoPriorityItem<DATA_T, PRIORITY_T>(data, priority, nextCount()));
    }

    /** topData
     *  Get the data value associated with top item on the queue.  This does
     *  not modify the queue.
     *
     *  @return DATA_T
     */
    DATA_T topData(void) {
        return getPriorityQueue().top().getData();
    }

    /** topPriority
     *  Get the priority value associated with top item on the queue.  This does
     *  not modify the queue.
     *
     *  @return PRIORITY_T
     */
    PRIORITY_T topPriority(void) {
        return getPriorityQueue().top().getPriority();
    }

    /** pop
     *  Remove the top item from the queue.
     */
    void pop(void) {
        getPriorityQueue().pop();
    }

    /** isEmpty
     *  Tests whether the queue is empty.
     *
     *  @return    bool
     */
    bool    isEmpty(void) {
        return getPriorityQueue().empty();
    }

  private:
    typedef std::priority_queue<FifoPriorityItem<DATA_T, PRIORITY_T>,
        std::vector<FifoPriorityItem<DATA_T, PRIORITY_T> >,
        FifoPriorityItemCompareFcn<DATA_T, PRIORITY_T> >     queue_T;
    
    uint64_T fCount;         /**< Counter that increments each time a push() happens. */
    queue_T fPriorityQueue;  /**< Priority queue from the standard library. */
    
    /** nextCount
     *  Returns the current count and then increments it.
     *
     *  @return    uint64_T count value.
     */
    uint64_T nextCount(void) {
        return fCount++;
    }

    /** getPriorityQueue
     *  Returns the contained priority queue (as a reference).
     *
     *  @return queue_T &
     */
    queue_T & getPriorityQueue(void) {
        return fPriorityQueue;
    }

    /** setPriorityQueue
     *  Sets the contained priority queue.
     *
     *  @param q   queue_T
     */
    void setPriorityQueue(queue_T q) {
        fPriorityQueue = q;
    }

    /** setCount
     *  Sets the internal counter.  Used by the constructor.
     *
     *  @param c   uint64_T
     */
    void setCount(uint64_T c) {
        fCount = c;
    }
};

#endif // FifoPriorityQueue_H
