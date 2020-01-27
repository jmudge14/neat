# neat

This is an implementation of the NEAT algorithm in Common Lisp, based directly on Kenneth Stanley's original (2002) paper, with a few tweaks (mostly around the reproductive processes, particularly in selecting which mutation stratgies to use, and some percentage of each generation is carried over into the next unconditionally.)

   http://nn.cs.utexas.edu/downloads/papers/stanley.ec02.pdf 

This is not a fork of the cl-neat package, but instead is implemented from scratch, primarily as an academic exercise.

Currently, this version converges on an xor solution in about 100 generations, which isn't particularly fast, but does show that it works. 

Future plans include the introduction of a "tournament"-style training routine, for self-reinforced learning. In this model, rather than applying a fitness test to individual genomes, they are compared by way of a tournament so that fitness is a relative number. (Comparison across generations is achieved by the aforementioned carryover percentage.) Ultimately, the intended goal is that this method could learn to play a board game against a human player without any input data other than a model of the board game itself.


### TODOs

* Code clean-up/refactor -- the code currently works, but it's messy.
* Implement tournament-style training mentioned above
* Bugfix: Xor sometimes gets "stuck", and instead of converging in about 100 generations it fails to converge within 1000 generations. Since this can be improved somewhat (or made significantly worse) by tweaking parameter weights, this may be a consequence of some difference between this and Stanley's original paper, which reports about 30 generations to converge.


### Author
Jack Mudge <jmudge14@gmail.com>_


## License

Public Domain

