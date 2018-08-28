/*
 * QueryManager.h
 *
 *  Created on: Jun 11, 2014
 *      Author: dkoes
 */

#ifndef QUERYMANAGER_H_
#define QUERYMANAGER_H_

#include <vector>
#include <string>
#include <boost/thread.hpp>
#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>
#include "MinimizationQuery.h"

using namespace boost;
using namespace std;

typedef boost::shared_ptr<MinimizationQuery> QueryPtr;

//an instance of this classes manages all the extant minimization queries
//each query is assigned a unique id for later reference
class QueryManager {
  private:
    unsigned nextID; //counter to generate unique IDs

    typedef unordered_map<unsigned, QueryPtr> QueryMap;
    QueryMap queries;

    boost::mutex mu;

    unsigned timeout; //seconds until purgeable

    MinimizationParameters minparm;
  public:

    QueryManager(unsigned numt, unsigned tout = 60 * 30)
        : nextID(1), timeout(tout) {
      minparm.nthreads = numt;
    }

    //add a query
    //first parse the text and return 0 if invalid
    //if oldqid is set, then deallocate/reuse it
    unsigned add(unsigned oldqid, stream_ptr io);

    QueryPtr get(unsigned qid);

    unsigned purgeOldQueries();

    void getCounts(unsigned& active, unsigned& inactive, unsigned& defunct);
    unsigned processedQueries() const {
      return nextID - 1;
    }

};

#endif /* QUERYMANAGER_H_ */
