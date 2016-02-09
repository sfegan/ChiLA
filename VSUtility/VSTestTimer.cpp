#include <cmath>
#include <iostream>

#include <VSTestTimer.hpp>

using namespace VERITAS;

#ifdef TEST_MAIN

double x;

int main(int argc, char** argv)
{
  VSTestTimer<VSTimerCoreIA32> a("a");
  //VSTestTimer<VSTimerCoreUNIX> a("a");
  a.start();
  sleep(5);
  a.stop();
  std::cerr << a << std::endl;

  VSTestTimer<VSTimerCoreUNIX> b("b");
  b.start();
  sleep(5);
  b.stop();
  std::cerr << b << std::endl;
}

#endif
