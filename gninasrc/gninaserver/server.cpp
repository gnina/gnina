//============================================================================
// Name        : server.cpp
// Author      : David Koes
// Created     : Jun 11, 2014
// Copyright   : 2014 University of Pittsburgh
// Description : Minimization server
//============================================================================

#include "server_common.h"
#include "CommandLine2/CommandLine.h"
#include "Logger.h"
#include "QueryManager.h"
#include "servercmds.h"

using namespace std;
using namespace boost;
using namespace boost::asio;
using namespace boost::asio::ip;

cl::opt<unsigned> port("port", cl::desc("port used by server"), cl::Required);
cl::opt<unsigned> maxConcurrent("max-concurrent-requests",
    cl::desc(
        "block further incoming requests after this amount of concurrency"),
    cl::init(16));
cl::opt<unsigned> minimizationThreads("threads",
    cl::desc("number of threads to use for minimization"),
    cl::init(max(1U, boost::thread::hardware_concurrency() / 2)));
cl::opt<string> logfile("logfile", cl::desc("file for logging information"));

typedef unordered_map<string, boost::shared_ptr<Command> > cmd_map;

static void process_request(stream_ptr s, cmd_map& cmap) {
  string cmd;
  //first line is command
  getline(*s, cmd);
  trim(cmd);

  try {
    if (cmap.count(cmd) > 0) {
      cmap[cmd]->execute(s);
    } else {
      *s << "ERROR\nInvalid command: " << cmd << "\n";
    }
  } catch (std::exception& e) {
    *s << "ERROR\nException " << e.what() << "\n";
  }
}

//periodically check for expired queries
static void thread_purge_old_queries(QueryManager *qmgr) {
  while (true) {
    boost::this_thread::sleep(posix_time::time_duration(0, 3, 0, 0));
    unsigned npurged = qmgr->purgeOldQueries();
  }
}

int main(int argc, char *argv[]) {
  cl::ParseCommandLineOptions(argc, argv);

  //setup log
  Logger log(logfile);
  QueryManager queries(minimizationThreads); //initialize query manager

  //command map
  cmd_map commands = assign::map_list_of("startmin",
      boost::shared_ptr<Command>(new StartMinimization(queries, log)))("cancel",
      boost::shared_ptr<Command>(new CancelMinimization(queries, log)))(
      "getscores", boost::shared_ptr<Command>(new GetScores(queries, log)))(
      "getjsonscores",
      boost::shared_ptr<Command>(new GetJSONScores(queries, log)))("getmol",
      boost::shared_ptr<Command>(new GetMol(queries, log)))("getmols",
      boost::shared_ptr<Command>(new GetMols(queries, log)))("getstatus",
      boost::shared_ptr<Command>(new GetStatus(queries, log)));

  //start listening
  io_service io_service;
  tcp::acceptor a(io_service, tcp::endpoint(tcp::v4(), port));

  cout << "Listening on port " << port << "\n";

  //start up cleanup thread
  boost::thread cleanup(thread_purge_old_queries, &queries);
  while (true) {
    stream_ptr s = stream_ptr(new tcp::iostream());
    a.accept(*s->rdbuf());
    boost::thread t(bind(process_request, s, commands));
  }

}
