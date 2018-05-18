Optimization {#cha-optimization}
============

While sampling from a posterior with Markov chains, BAT will store the
maximum value of the posterior that it has seen (and the corresponing
point in parameter space). But the sampling algorithm is designed to
sample, not to find a maximum value. To find the maximum point, BAT
includes optimization methods. These are called via
`BCIntegrate::FindMode`, which takes optional arguments: You may give
it a starting point; if you do not specify one, BAT will start from
the current maximum point. You can also specify the optimization
method when you call the function; or you can specify it in advance
via `BCIntegrate::SetOptimzationMethod`.

By default, BAT replaces its currently held maximum point only if the
newly found one has a greater posterior. You can tell BAT to replace
the currently held one regardless of whether optimization improves
upon the currently stored one via the function
`BCIntegrate::SetFlagIgnorePrevOptimization`.

# Minuit {#sec-minuit}

BAT uses Minuit through ROOT's interface to it. To use it, set the
optimization method to `BCIntegrate::kOptMinuit`.

Minuit is what is known as a gradient follower---it moves from a
starting point in the direction that increases the posterior until it
finds a maximum. Minuit is much better at finding the exact location
of a maximum than sampling with Markov chains is, but it does not
gaurantee that this maximum is the global maximum or only a local
one. Sampling with Markov chains can better identify the region of the
global maximum than Minuit can. So we recommend that you first
marginalize your model and then call `FindMode`.

# Simulated Annealing {#sec-simulated-annealing}

BAT provides a simulated annealing algorithm for optimization. To use
it, set the optimization method to `BCIntegrate::kOptSimAnn`.

This algorithm is similar to the Metropolis-Hastings sampling one in
that involves proposal of new points to move to randomly in a
neighborhood of the current point. But the neighborhood and the
acceptance criteria are regulated in such a way as to encourage motion
towards the global maximum.

Using BAT, simulated annealing will in general not find the maximum
point as precisely or rapidly as Minuit, but it can more reliably find
the global maximum instead of a local maximum. We therefore recommend
calling `FindMode` with `BCIntegrate::kOptMinuit` after calling it
with `kOptSimAnn`. This is similar to first sampling and then
optimizing with Minuit; but simulated annealing is less
calculationally expensive than sampling---that is, it will call your
likelihood less often. However, simulated annealing is not a sampling
technique and will not provide output samples for further offline use.

There are several parameters that govern the running of simulated
annealing. Please consult the code documentation and C. Brachem's
_Implementation and test of a simulated annealing algorithm in the
Bayesian Analysis Tookit_ (2009) for how to set them. The most
important option you may change is the proposal function. This is set
via `BCIntegrate::SetSASchedule`. The options are
`BCIntegrate::kSACauchy`, `BCIntegrate::kSABoltzmann`, and
`BCIntegrate::kSACustom`. (The latter requires you to set your own
proposal function.)
