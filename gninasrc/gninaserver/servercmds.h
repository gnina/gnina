/*
 * servercmds.h
 *
 *  Created on: Jun 11, 2014
 *      Author: dkoes
 *
 *  Classes that process various server commands.
 */

#ifndef SERVERCMDS_H_
#define SERVERCMDS_H_

#include "Logger.h"
#include "QueryManager.h"
#include "server_common.h"
#include <fstream>

class Command {
  protected:

    Logger& log;

  public:
    Command(Logger& l)
        : log(l) {
    }
    virtual ~Command() {
    }

    virtual void execute(stream_ptr io) = 0;

};

//start a minimization which gets added to the query manager
class StartMinimization : public Command {
    QueryManager& qmgr;

  public:
    StartMinimization(QueryManager& q, Logger& l)
        : Command(l), qmgr(q) {
    }

    void execute(stream_ptr io) {
      //next line is an old qid
      string str;
      getline(*io, str);
      trim(str);
      unsigned oldqid = atoi(str.c_str());

      unsigned qid = qmgr.add(oldqid, io);
      //return the query id (zero if there's a problem)

      log.log("startmin %d %d\n", oldqid, qid);
      *io << qid << "\n";
      io->flush();
      //do not close io, need to finish reading
    }
};

//cancel a minimization
class CancelMinimization : public Command {
    QueryManager& qmgr;

  public:
    CancelMinimization(QueryManager& q, Logger& l)
        : Command(l), qmgr(q) {
    }

    void execute(stream_ptr io) {
      //next line is an old qid
      string str;
      getline(*io, str);
      trim(str);
      unsigned oldqid = atoi(str.c_str());
      log.log("cancel %d\n", oldqid);
      QueryPtr query = qmgr.get(oldqid);
      if (query) {
        query->cancel();
      }
      io->close();
    }
};

//get the text scoring output of the specified minimization
class GetScores : public Command {
    QueryManager& qmgr;

  public:
    GetScores(QueryManager& q, Logger& l)
        : Command(l), qmgr(q) {
    }

    void execute(stream_ptr io) {
      //query id followed by filter params
      MinimizationFilters filters;
      unsigned qid = 0;
      *io >> qid;
      filters.read(*io);
      QueryPtr query = qmgr.get(qid);
      if (query) {
        query->outputData(filters, *io);
      }
      io->close();
    }
};

//get the json scoring output of the specified minimization
class GetJSONScores : public Command {
    QueryManager& qmgr;

  public:
    GetJSONScores(QueryManager& q, Logger& l)
        : Command(l), qmgr(q) {
    }

    void execute(stream_ptr io) {
      //query id followed by filter params
      MinimizationFilters filters;
      unsigned qid = 0, draw = 0;
      *io >> qid;
      *io >> draw; //datatable draw code
      filters.read(*io);
      QueryPtr query = qmgr.get(qid);
      if (query) {
        query->outputJSONData(filters, draw, *io);
      }
      io->close();
    }
};

//return a single requested minimized molecule structure
class GetMol : public Command {
    QueryManager& qmgr;

  public:
    GetMol(QueryManager& q, Logger& l)
        : Command(l), qmgr(q) {
    }

    void execute(stream_ptr io) {
      //first line query id and mol id
      unsigned qid = 0, molid = 0;
      *io >> qid;
      *io >> molid;
      QueryPtr query = qmgr.get(qid);
      if (query) {
        query->outputMol(molid, *io);
      }
      io->close();
    }
};

//return all the minimized molecule structures - gzipped
class GetMols : public Command {
    QueryManager& qmgr;

  public:
    GetMols(QueryManager& q, Logger& l)
        : Command(l), qmgr(q) {
    }

    void execute(stream_ptr io) {
      //query id followed by filter params
      MinimizationFilters filters;
      unsigned qid = 0;
      *io >> qid;
      filters.read(*io);
      QueryPtr query = qmgr.get(qid);
      if (query) {
        query->outputMols(filters, *io);
      }
      io->close();
    }
};

//return current status
class GetStatus : public Command {
    QueryManager& qmgr;

  public:
    GetStatus(QueryManager& q, Logger& l)
        : Command(l), qmgr(q) {
    }

    void execute(stream_ptr io) {
      unsigned active, inactive, defunct;
      qmgr.getCounts(active, inactive, defunct);
      double load = 0;

      ifstream ldfile("/proc/loadavg");
      ldfile >> load;

      *io << "Active " << active << "\nInactive " << inactive << "\nDefunct "
          << defunct << "\nLoad " << load << "\n";
      io->close();
    }
};

#endif /* SERVERCMDS_H_ */
