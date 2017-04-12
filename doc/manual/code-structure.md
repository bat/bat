Structure of the Code {#cha-code-structure}
=====================

@todo 1-2 pages, inheritance

[TOC]

BAT is object-oriented and uses inheritance for code reuse. The single
most important hierarchy is that of `BCModel`. In most applications, a
user would implement a model `MyModel` by either inheriting from
`BCModel` directly or creating an instance of one of the @ref
cha-predefined-models and thus would inherit from `BCModel`
indirectly.

@dot
digraph G {
    node [shape=box];
    a [ label="BCEngineMCMC" URL="@ref BCEngineMCMC" ];
    b [ label="BCIntegrate"  URL="@ref BCIntegrate" ];
    c [ label="BCModel"      URL="@ref BCModel" ];
    a -> b; b -> c;
}
@enddot

@htmlonly
The boxes on the diagram are links to the reference guide for the respective classes.
@endhtmlonly

This structure has the benefit that `MyModel` has native access to all
things related to Markov chains due to `BCEngineMCMC` and one can
immediately optimize the posterior with the methods in
`BCIntegrate`. However, one has to keep in mind that multiple Markov
chains for the *same* instance of `MyModel` share the entire state of
the object. This can lead to complications when
`MyModel::LogLikelihood` is called concurrently; see @ref
cha-threading for more details.

@section sec-code-variables Variables

@htmlonly
<br/><br/>
Return to the <a href="cha-predefined-models.html">previous</a> section or go to the
<a href="cha-output.html">next</a> section.
@endhtmlonly
