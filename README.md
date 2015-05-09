# hcc_dclinde

A tool to infer discrete High-Order Dynamic Bayesian Network (HO-DBN)
with hidden common cause(s) for gene regulatory network from discrete
time series data, where the expression of a small but unknown number
of nodes is hidden.

====================================================================
Introduction to hcc_dclinde
====================================================================

hcc_dclinde is an algorithm to recover the causal GRN with delays in
the regulatory links from time series expression data, where a small
but unknown number of nodes are hidden, i.e. without expression data.

We assume that there is a d-th order (the maximum delay is d)
stationary HO-DBN that is of interest, where a small but unknown
number of common cause(s) are not observed. Each hidden variable has
only observed variables as children and parents, with at least two
children and possibly no parents. We also make the simplifying
assumption that the children of unobserved variable(s) are not linked
to each other, because it is difficult to differentiate whether the
high association between two children are solely due to the hidden
common cause, or due to both the hidden common cause and direct link
between the two children. As the prior network is difficult to learn,
and the transition network is of the main interest, our objective is
to infer the transition network of the HO-DBN from the discrete time
series of the observed variables. Since it is difficult to obtain very
long time series, hcc_dclinde is also capable of utilizing multiple
short time series, which are not necessarily of the same length
(e.g. obtained from replicate experiments). The programs have been
tested in Debian and Ubuntu Linux distributions, and should build and
work with minimal tweaking in other Linux distributions.

==== Building gsl

Before buliding hcc_dclinde, you need to build gsl, which is the GNU
Scientific Library, which contains routines used by hcc_dclinde.

gsl-1.16.tar.gz is included for convenience, and we have tested using
gsl-1.16. First extract it by:

    tar xzvf gsl-1.16.tar.gz

Then to build it, typically you can use:

    cd gsl-1.16/
    ./configure
    make
    make install

But refer to gsl-1.16/README for trouble shooting.

==== Building hcc_dclinde

Having built gsl, you could build hcc_dclinde by using the provided makefile:

    make

You may also build it directly, using:

    g++ -O2 -Wall -c -o globalMIT.o globalmit_source/globalMIT.cpp
    g++ -O2 -Wall -c -o globalMIT_fixedOrder.o globalmit_source/globalMIT_fixedOrder.cpp
    gcc -Wall -O2 -static hcc_dclinde.c parse_option.c dtsv.c globalMIT.o globalMIT_fixedOrder.o mt19937ar_lib.o -lgsl -lgslcblas -lm -lstdc++

==== Usage:

The usage of hcc_dclinde is:

     Usage: ./hcc_dclinde [-?] -data item1 [item2 ...] [-seed str] [-st1 real]
     [-st2 real] [-min.delay int] [-max.delay int] [-max.n int] [-pruning str]
     [-keep.zero.delay] [-one.delay] [-no.dup] [-istart int] [-iend int]
     [-alpha real] [-no_self_link] [-no.cpt] [-out.grn str] [-out.init.grn str]
     [-init.only] [-method1 str] [-method2 str] [-eV real] [-eV.tolerance real]
     [-eV.threshold real] [-nh.states int] [-n.EM.iter int]
     [-n.no.change.restart int] [-mp int] [-hp.st real]
     
     Description of the options:
       -?:  Showing the usage.
     
       -data [REQUIRED]:  File name(s) of the input expression data (tab/space
         separated) small non-negative integers for discrete data, each row is a
         time point, each column is a gene. Each file should have the same number
         of columns, but may have different number of rows.
     
       -seed:  The seed of the pseudo random number generator for EM initial
         distribution. If contains '-' (.e.g negative number), use current time.
         Default "0"
     
       -st1:  Score threshold for stage 1 of D-CLINDE. Score is -log10(p), where p
         is the intended p-value. Default 2.
     
       -st2:  Score threshold for stage 2 of D-CLINDE. Similar to st1. Default 2.
     
       -min.delay:  Minimum delay (lags) to use in the inference (not applicable to
         GlobalMIT+/*), default 1.
     
       -max.delay:  Maximum delay (lags) to use in the inference, default 4.
     
       -max.n:  Maximum number of parents to condition on in stage 2 of D-CLINDE,
         default 2.
     
       -pruning:  Pruning strategy in stage 2 of D-CLINDE, can be all (consider all
         neighbors of the two vertices when doing conditional test of the link
         between the two) or common (consider only common neighbors of the two
         vertices). Default "all"
     
       -keep.zero.delay:  To keep links with zero delay after stage 2 of D-CLINDE
         and before remove dups. Default false.
     
       -one.delay:  To keep only one delay (the one with best score, smallest
         delay) for each link after stage 1 of D-CLINDE. Default false.
     
       -no.dup:  Remove duplicate links. To keep only one delay (the one with best
         score, smallest delay) for each link after stage 2 of D-CLINDE. Default
         false.
     
       -istart:  The starting 0-based index of gene for which to find parents. Only
         for GlobalMIT+/* in initial GRN, default 0.
     
       -iend:  The ending (inclusive) 0-based index of gene for which to find
         parents. Only for GlobalMIT+/* in initial GRN, default -1, which means
         n_genes-1.
     
       -alpha:  Significance level for GlobalMIT+/*. Default 0.999.
     
       -no_self_link:  To disallow self link (with delay) for GlobalMIT+/*, default
         false.
     
       -no.cpt:  Do not print the conditional distributions for the final GRN.
         Default false.
     
       -out.grn:  Name of file to output the final GRN. Default NULL, which means
         do not output to file.
     
       -out.init.grn:  Name of file to output the init GRN. Default NULL, which
         means do not output to file.
     
       -init.only:  To learn only the initial GRN, which does not consider possible
         hidden common causes. Default false.
     
       -method1:  Method for learning initial GRN, can be either "dclinde" for
         D-CLINDE, "mit" for GlobalMIT+, or "mit_fast" for GlobalMIT*. Default
         "mit_fast"
     
       -method2:  Method for learning sub-GRN after potential hidden common
         cause(s) are identified and estimated, can be either "dclinde" for
         D-CLINDE, "mit" for GlobalMIT+, or "mit_fast" for GlobalMIT*. Default
         "mit_fast"
     
       -eV:  Expected Bias if > 0. If <= 0, the expected bias will be estimated to
         be the median of the bias of the initial GRN. If bias of a gene is > eV +
         eV.tolerance or < eV - eV.tolerance, that gene will be considered for
         having hidden common node. Default 0.65.
     
       -eV.tolerance:  If bias of a gene is > eV + eV.tolerance or < eV -
         eV.tolerance, that gene will be considered for having hidden common node.
         Default 0.05.
     
       -eV.threshold:  Score threshold of association for clustering using G2 test.
         Association score is -log_10(p-value). Similar to st1 and st2. Default
         2.3.
     
       -nh.states:  Number of states for hidden cause(s). If < 2, use the maximum
         number of states of children. Default 0.
     
       -n.EM.iter:  Number of iterations for EM in estimating hidden cause(s).
         Default 100.
     
       -n.no.change.restart:  In EM, if the most probable hidden states do not
         change after n.no.change.restart iterations, restart EM. Default 3.
     
       -mp:  The maximum number of potential parents to find for hidden cause after
         the first EM. Default 3.
     
       -hp.st:  Score threshold of association using G2 test for finding potential
         parent for hidden cause. Association score is -log_10(p-value). Similar to
         st1 and st2. Default 2.
         

==== Example usage:

An example use of hcc_dclinde is for n=50 sample data using D-CLINDE
for initial GRN and re-learning final GRN:

    ./hcc_dclinde -eV 0 -method1 dclinde -method2 dclinde -data sample_data/r2_s*.txt -out.grn test_grn.txt > output.txt

Now output.txt contains some messages and the GRN after various stages
of hcc_dclinde, and test_grn.txt contains the final GRN in format
acceptable for the comparison tool grn_cmp_hcc.

Note that for n=50, the default method mit_fast may still be too slow
on a slow computer.

====================================================================
Comparison of two GRNs (if applicable)
====================================================================

If the true GRN is known, you may assess how close the predicted GRN
is to the true GRN by using grn_cmp_hcc.

==== Building grn_cmp_hcc:

You may build it by using the provided makefile:

    make grn_cmp_hcc

Or you may build it directly by:

    gcc -Wall -O3 grn_cmp_hcc.c parse_option.c -o grn_cmp_hcc -lm

==== Usage of grn_cmp_hcc:

     Usage: ./grn_cmp_hcc [-?] -p str -t str [-n int] [-no.effect] [-no.delay] [-v]

     Description of the options:
       -?:  Showing the usage.
     
       -p [REQUIRED]:  File name of the predicted GRN. Each line in the file
         represents an edge, which consists of a 0-based 'from' index, a 0-based
         'to' index, a delay, and the effect, separated by space..
     
       -t [REQUIRED]:  File name of the true GRN. Each line in the file represents
         an edge, which consists of a 0-based 'from' index, a 0-based 'to' index, a
         delay, and the effect, separated by space..
     
       -n:  Number of non-hidden genes, unspecified if <= 0. If an index is >=
         this, it is regarded as a hidden node. Default 0.
     
       -no.effect:  No effect.
     
       -no.delay:  No delay.
     
       -v:  Verbose mode.
     

The parameter -n specifies the number of observed genes, which is the
number of columns in the sample data. -no.effect should be used for
discrete GRN, as in this case.


==== Example usage:

An example use of grn_cmp_hcc is:

    ./grn_cmp_hcc -p test_grn.txt -t sample_data/r2_grn.txt -n 50 -no.effect

====================================================================
R Scripts for Generation of Synthetic Data (Optional)
====================================================================

The following functions in synthetic.R can help generate three
different types of synthetic data for testing purposes. Of course R
need to be installed first.

Useful functions:

synthetic_small_cases_e.R:

    "gen.small.hidden.cases()": small case with one hidden node

    "gen.small.non.hidden.cases()": small case with no hidden node

    "gen.misc.cases()": larger case with more than one hidden node

Be warned that a many synthetic data files will be generated, and it
may take some time because R is not very fast, and the script are not
written in a particularly efficient way.

You may tweak the functions to generate cases with other parameters.

====================================================================

====================================================================
Licence
====================================================================
The following files are written by Peter Lo (me, email:
peter19852001@yahoo.com.hk), and I hereby grant the permission of
anyone to use these files for any non-commercial purposes, with not
warranty.:

hcc_dclinde.c
parse_option.{c,h}
grn_cmp_hcc.c
dtsv.{c,h}

And the files in globalmit_source are originally written by Vinh
Nguyen, and modified by Peter Lo to use in hcc_dclinde, with approvabl
(by email) from the original author. You should seek approval from the
original author to modify the code.

====================================================================