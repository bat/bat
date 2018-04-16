How to run unit tests
=====================

1. Run all tests: in $(top_srcdir), `make check`. This builds BAT if necessary.
2. Run individual test: in $(top_srcdir), `make check TESTS=BCModel_TEST`
3. Run a test in a debugger: in $(top_srcdir), `make check TESTS=` to build all tests. In `test/`, there are only scripts that setup the environment. The actual executable is in `test/.libs`. To debug it with for example `eclipse`, choose `test/.libs/lt-BCModel_TEST`.

How to add a new unit test
==========================

1. Open `Makefile.am`, add `foo.TEST` to `TESTS` and add the sources as `foo_TEST_SOURCES = foo_TEST.cxx`.
2. Create `foo_TEST.cxx`, enter this skeleton:

```cpp
/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include <test.h>

using namespace test;

class FooTest :
    public TestCase
{
public:
    FooTest() :
        TestCase("Foo test")
    {
    }

    virtual void run() const
    {
        TEST_CHECK_EQUAL(1+1, 3);
    }
} fooTest;
```

Further test macros are defined in `test.h`.

3. Run the test and make sure it *fails* the first time so you are sure it is being run.
