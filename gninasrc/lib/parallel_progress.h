/*

 Copyright (c) 2006-2010, The Scripps Research Institute

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 Author: Dr. Oleg Trott <ot14@columbia.edu>, 
 The Olson Lab, 
 The Scripps Research Institute

 */

#ifndef VINA_PARALLEL_PROGRESS_H
#define VINA_PARALLEL_PROGRESS_H

#include <boost/thread/mutex.hpp>
#include <boost/timer/timer.hpp>
#include <boost/noncopyable.hpp>
#include <boost/cstdint.hpp>  // for uintmax_t
#include <iostream>           // for ostream, cout, etc
#include <string>

#include "incrementable.h"

//taken from boost/progress.hpp as it has been deprecated
class progress_display : private boost::noncopyable
{
 public:
  explicit progress_display( unsigned long expected_count_,
                             std::ostream & os = std::cout,
                             const std::string & s1 = "\n", //leading strings
                             const std::string & s2 = "",
                             const std::string & s3 = "" )
   // os is hint; implementation may ignore, particularly in embedded systems
   : noncopyable(), m_os(os), m_s1(s1), m_s2(s2), m_s3(s3) { restart(expected_count_); }

  void           restart( unsigned long expected_count_ )
  //  Effects: display appropriate scale
  //  Postconditions: count()==0, expected_count()==expected_count_
  {
    _count = _next_tic_count = _tic = 0;
    _expected_count = expected_count_;

    m_os << m_s1 << "0%   10   20   30   40   50   60   70   80   90   100%\n"
         << m_s2 << "|----|----|----|----|----|----|----|----|----|----|"
         << std::endl  // endl implies flush, which ensures display
         << m_s3;
    if ( !_expected_count ) _expected_count = 1;  // prevent divide by zero
  } // restart

  unsigned long  operator+=( unsigned long increment )
  //  Effects: Display appropriate progress tic if needed.
  //  Postconditions: count()== original count() + increment
  //  Returns: count().
  {
    if ( (_count += increment) >= _next_tic_count ) { display_tic(); }
    return _count;
  }

  unsigned long  operator++()           { return operator+=( 1 ); }
  unsigned long  count() const          { return _count; }
  unsigned long  expected_count() const { return _expected_count; }

  private:
  std::ostream &     m_os;  // may not be present in all imps
  const std::string  m_s1;  // string is more general, safer than
  const std::string  m_s2;  //  const char *, and efficiency or size are
  const std::string  m_s3;  //  not issues

  unsigned long _count, _expected_count, _next_tic_count;
  unsigned int  _tic;
  void display_tic()
  {
    // use of floating point ensures that both large and small counts
    // work correctly.  static_cast<>() is also used several places
    // to suppress spurious compiler warnings.
    unsigned int tics_needed = static_cast<unsigned int>((static_cast<double>(_count)
        / static_cast<double>(_expected_count)) * 50.0);
    do { m_os << '*' << std::flush; } while ( ++_tic < tics_needed );
    _next_tic_count =
      static_cast<unsigned long>((_tic/50.0) * static_cast<double>(_expected_count));
    if ( _count == _expected_count ) {
      if ( _tic < 51 ) m_os << '*';
      m_os << std::endl;
      }
  } // display_tic
};

struct parallel_progress : public incrementable {
    parallel_progress()
        : p(NULL) {
    }
    void init(unsigned long n) {
      p = new progress_display(n);
    }
    void operator++() {
      if (p) {
        boost::mutex::scoped_lock self_lk(self);
        ++(*p);
      }
    }
    virtual ~parallel_progress() {
      delete p;
    }
  private:
    boost::mutex self;
    progress_display* p;
};

#endif

