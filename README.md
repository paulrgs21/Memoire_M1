SECTION 4
In the sensory analysis context, there is generally no censoring in the experiments
since the “panelists” stop by themselves the tasting procedure. That is to say, each
trajectory Si is observed until it reaches an absorbing state -> nous qu'est ce qu'on considère étant donné nos data ?

SECTION 5
We suppose, as in Lecuelle et al. (2018) and Cardot et al. (2019), that the distribu-
tions of the sojourn times only depend on the current state, meaning that there is no
anticipation of the next state of the chain ->  pas sûr pour trajectoires professionnelles !

doParallel libraries.

The choice of the sample size and of
the number D of states is chosen to mimick the usual empirical framework in sensory
analysis studies -> nous qu'est ce qu'on considère étant donné nos data 

There is no explicit solution for the maximum likelihood estimators and we obtain
the estimated values of the parameters with the numerical Nelder-Mead optimization
routines in R -> à noter pour MLE pour gamma et weibull

x comprendre les ddl fin page 14

A precise description of the algorithm used
to generate transition matrices is given in the supplementary file -> ça se fait sans ?

The values of the parameters characterizing the Gamma and the Weibull distribution
are the values estimated on the dataset of chocolates (see Sect. 6). -> EMV OU ALGO NELDER MEAD

In this dataset, the number of states D was equal to 10. We also simulated trajectories of SMP with D = 4
and D = 7 states by considering only the simulated transition probabilities and the
sojourn time distributions related to the last four and last seven states. As far as the
absorbing case is concerned, we consider D A = D + 1 attributes. -> PAS COMPRIS

We present in the next subsections the results considering a Gamma distribution
for the sojourn times. The results for the Weibull distribution, which lead to similar
conclusions, are presented in the supplementary file. -> encore ce file, on s'en fout ?



5.2 POWER OF THE TEST -> balek ? en vrai c'est important car "L'efficacité d'un test statistique est généralement évaluée selon deux critères :
    -Le respect du niveau nominal (γ=0.05γ=0.05 ou 0.10.1) : un bon test ne doit pas rejeter H0H0​ plus souvent que prévu sous l'hypothèse nulle.
    -La puissance : un bon test doit détecter une différence réelle lorsque celle-ci existe (c'est-à-dire rejeter H0H0​ quand H1H1​ est vraie)."

ω ε2 = (1 − ε2)ω0 + ε2ω1 ?? WTF ? ce sont les paramètres estimés et pas le nombre de paramètres ? 

As the matrices P1 and P2 are obtained by simulations as explained in the previous
Section and the parameters of the Gamma and the Weibull law are obtained thanks to
the dataset of chocolates, it is not easy to quantify how the parameters are different. -> pas plutôt P0 et P1 ? comment on fait pour estimer les paramètres si on n'a pas les données ?

selon le milieu de la page 18, les points 2 et 3 correspondent aux tests partiels donc on peut les avoid


Si on observe comme eux sur données simulées que permutations est le meilleur test, et que khi-2 est fine avec n grand, D/n petit (et donc calculs très rapides), pas besoin de faire les 3 méthodes sur les données réelles non ?


SECTION 6

à faire sur nos données réelles : comparaison entre weibull, gamma, exp pour savoir quelle distribution suivent réellement nos temps de séjour (voir table 6)
