Multithreading and Thread Safety {#cha-threading}
================================

In a demanding analysis, a single evaluation of the posterior may take seconds
or more. In an MCMC analysis, this is then the bottle neck. The overall time to
solution can then be reduced by running multiple chains on multiple threads.

Assuming BAT is configured with parallelization enabled (see the @ref cha-install), the number of threads can be selected at runtime
without recompilation with `./program OMP_NUM_THREADS=N`. Due to the
overhead from thread creation, the sampling is faster only if the likelihood is
sufficiently slow to compute. As a rule of thumb, if the unparallelized sampling
takes minutes or even hours, the parallelized version should get close to the
maximum speedup given by the number of cores. For very simple likelihoods (say one millisecond per evaluation), the
overhead from synchronizing threads usually slows down the entire process. We invite the user to experiment which number of threads gives the best performance.

If BAT is compiled with support for threads and the number of threads
is not set explicitly, it is implementation dependent whether only one
thread is created or as many as there are available cores. To enforce serial execution, simply set `OMP_NUM_THREADS=1`.

For guidance, tests showed that it is beneficial to have as many threads as your
computer provides. Working on a quad core machine with hyper-threading, we
observed a speedup factor of 3.5 -- 3.9 with 8 and 16 chains on 8 cores for a likelihood that takes more than a second to evaluate.

@warning The results of the sampling become nonsense if multiple threads
are used but the likelihood itself is not thread safe.

A simple example of a non-thread-safe likelihood is
@code{.cpp}
double MyModel::LogLikelihood(const std::vector<double>& parameters)
{
  // assign member variable
  this->member = parameters[0];

  // value of member used in MyModel::Method
  return this->Method();
}

double MyModel::Method()
{
  return -member * member;
}
@endcode

If the method `MyModel::Method` is called simultaneously from different threads
with different `parameters`, its return value depends on the arbitrary call
order of individual threads. There is a race condition on `member` because it is
read and written simultaneously by multiple threads so the result of
`LogLikelihood` is not deterministic.
There are several solutions.

## Avoid state

In this contrived example, it is easiest to avoid using any member variable or method call and do everything inside `LogLikelihood`
@code{.cpp}
double MyModel::LogLikelihood(const std::vector<double>& parameters)
{
  return -parameter[0] * parameter[0];
}
@endcode

## Independent copies of state

In practice, it may be impractical or even inevitable to maintain state and call into other
functions that require state. A clean solution is to move `Method` into
`OtherClass`, and to keep an independent copy of `OtherClass` for each thread or
each chain. For this purpose, we provide the hook
`BCEngineMCMC::MCMCUserInitialize` that is called before any chain attempts to
evaluate the likelihood. For example:

@code{.cpp}
class SomeState
{
public:
  double Method() { ... }
private:
  double t;
}

class MyModel
{
public:
  ...
  virtual void MCMCUserInitialize()
  {
    state.resize(GetNChains());
  }

  virtual double LogLikelihood(const std::vector<double>& parameters)
  {
    SomeState& s = state.at(GetCurrentChain());
    s.t = parameters[0];
    return s.Method();
  }

private:
  std::vector<SomeState> state;
};
@endcode

@note In BAT, only the Markov chain sampling uses multiple threads. During optimization or any algorithm, only one thread is active. In the above example, `BCEngineMCMC::GetCurrentChain` returns 0 in such a context.
