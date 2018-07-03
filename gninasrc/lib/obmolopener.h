/*
 * obmolopener.h
 *
 *  Created on: Jan 10, 2013
 *      Author: dkoes
 *
 * This is a wrapper class for cleanly opening possibly gzipped molecular
 * data files.  It will delete and invalidate any streams that are still
 * open when it goes out of scope.
 */

#ifndef OBMOLOPENER_H_
#define OBMOLOPENER_H_

#include <openbabel/obconversion.h>
#include <fstream>
#include <vector>

class obmol_opener {
  public:

    obmol_opener() {
    }
    virtual ~obmol_opener();

    void openForInput(OpenBabel::OBConversion& conv, const std::string& name);
    void openForOutput(OpenBabel::OBConversion& conv, const std::string& name);

    void clear();

  private:

    obmol_opener& operator=(const obmol_opener& rhs); //you don't want to do this
    std::vector<std::ios*> streams;
};

#endif /* OBMOLOPENER_H_ */
