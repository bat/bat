/*
 * Copyright (C) 2012,  Danny van Dyk and Frederik Beaujean
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef __BC_TEST__H
#define __BC_TEST__H 1

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>

namespace test
{
    namespace implementation
    {
        template <typename T_>
        struct DoStringify
        {
            static std::string stringify(const T_ & x, unsigned precision)
            {
                std::stringstream ss;
                ss.precision(precision);
                ss << x;

                return ss.str();
            }
        };

        template <>
        struct DoStringify<std::string>
        {
            static std::string stringify(const std::string & x, unsigned)
            {
                return x;
            }
        };
    }

    template <typename T_>
    std::string stringify(const T_ & x, unsigned precision = 10)
    {
        return implementation::DoStringify<T_>::stringify(x, precision);
    }

    /*!
     * Print a range of a STL container or ordinary numerical array into a string in parentheses
     */
    template <typename Iterator_>
    std::string print_vector(const Iterator_ & begin, const Iterator_ & end, unsigned precision = 10)
    {
        std::stringstream ss;
        ss.precision(precision);
        ss << '(';

        for (Iterator_ i = begin ; i != end ; ++i)
        {
            ss << ' ' << *i;
        }

        ss << " )";

        return ss.str();
    }

    /*!
     * Print a STL container into a string.
     * Elements are surrounded by parentheses.
     * @param container The STL container from which elements are read.
     * @param precision String output precision.
     * @return String representation.
     */
    template <typename Container_>
    std::string print_container(const Container_ & container, unsigned precision = 10)
    {
        return print_vector(container.begin(), container.end(), precision);
    }

    class TestCase
    {
        private:
            std::string _name;

        public:
            TestCase(const std::string & name);

            virtual ~TestCase();

            std::string name() const;

            virtual void run() const = 0;
    };

    class TestCaseFailedException
    {
        private:
            int _line;

            std::string _file;

            std::string _reason;

        public:
            TestCaseFailedException(int line, const std::string & file, const std::string & reason);

            const std::string & reason() const;

            std::string where() const;
    };

#define TEST_SECTION(name, body) \
    do \
    { \
        std::cout << name << "> begins" << std::endl; \
        body \
        std::cout << name << "> ends" << std::endl; \
    } \
    while (false)

#define TEST_CHECK(a) \
    do \
    { \
        if (! (a)) \
            throw TestCaseFailedException(__LINE__, __FILE__, "'" #a "' is false"); \
    } \
    while (false)

#define TEST_CHECK_EQUAL(a, b) \
    do \
    { \
        if (! ((a) == (b))) \
            throw TestCaseFailedException(__LINE__, __FILE__, "'" #a "' is not equal to '" #b "'"); \
    } \
    while (false)

#define TEST_CHECK_FAILED(s) \
    do \
    { \
        throw TestCaseFailedException(__LINE__, __FILE__, s); \
    } \
    while (false)

#define TEST_CHECK_NEARLY_EQUAL(a, b, eps) \
    do \
    { \
        double a_val = (a); \
        double b_val = (b); \
        if (std::abs((a_val - b_val)) <= eps) \
            break; \
        else \
            throw TestCaseFailedException(__LINE__, __FILE__, \
                    "'" #a "' = " + stringify(a_val, 16) + " is not nearly-equal to '" #b "' = " + stringify(b_val, 16) + " within '" + stringify(eps) + "'" \
                    + ", difference is '" + stringify(a_val - b_val) + "'"); \
    } \
    while (false)

#define TEST_CHECK_RELATIVE_ERROR(a, b, eps) \
    do \
    { \
        double a_val = (a); \
        double b_val = (b); \
        if (std::sqrt(std::fabs(a_val)) < std::numeric_limits<double>::epsilon()) \
            throw TestCaseFailedException(__LINE__, __FILE__, \
                    "'" #a "' has been evaluated to the zero within computational accuracy, result = " + stringify(a_val)); \
         \
        if (std::sqrt(std::fabs(b_val)) < std::numeric_limits<double>::epsilon()) \
            throw TestCaseFailedException(__LINE__, __FILE__, \
                    "'" #b "' has been evaluated to the zero within computational accuracy, result = " + stringify(b_val)); \
         \
        if (((std::abs((a_val - b_val) / a_val)) <= eps) && ((std::abs((a_val - b_val) / b_val)) <= eps)) \
            break; \
        else \
            throw TestCaseFailedException(__LINE__, __FILE__, \
                    "One relative error of '" #a "' = '" + stringify(a_val, 16) + "' and '" #b "' = '" + stringify(b_val, 16) + "' is greater than " + stringify(eps) + ". The results are " + \
                    stringify(std::abs((a_val - b_val) / a_val)) + " and " + stringify(std::abs((a_val - b_val) / b_val))); \
    } \
    while (false)

#define TEST_CHECK_NO_THROW(expression) \
    do \
    { \
        try \
        { \
            expression; \
        } \
        catch (...) \
        { \
            throw TestCaseFailedException(__LINE__, __FILE__, \
                    "Caught unexpected exception in '" #expression "'"); \
        } \
    } \
    while (false)

#define TEST_CHECK_THROWS(exception, expression) \
    do \
    { \
        try \
        { \
            expression; \
        } \
        catch (exception & e) \
        { \
            break; \
        } \
        catch (...) \
        { \
            throw TestCaseFailedException(__LINE__, __FILE__, \
                    "Caught unexpected exception when expecting " #exception " in '" #expression "'"); \
        } \
        \
        throw TestCaseFailedException(__LINE__, __FILE__, \
                "Caught no exception in " #expression " when expecting '" #exception "'"); \
    } \
    while (false)
}

#endif
