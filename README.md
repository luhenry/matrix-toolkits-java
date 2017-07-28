[![Build Status](https://travis-ci.org/andreas-solti/matrix-toolkits-java.svg?branch=master)](https://travis-ci.org/andreas-solti/matrix-toolkits-java)
[![Coverage Status](https://coveralls.io/repos/fommil/matrix-toolkits-java/badge.svg?branch=master)](https://coveralls.io/r/fommil/matrix-toolkits-java?branch=master)
[![Maven Central](https://maven-badges.herokuapp.com/maven-central/io.github.andreas-solti.matrix-toolkits-java/mtj/badge.svg)](https://maven-badges.herokuapp.com/maven-central/io.github.andreas-solti.matrix-toolkits-java/mtj)
[![Javadoc](https://javadoc-emblem.rhcloud.com/doc/io.github.andreas-solti.matrix-toolkits-java/mtj/badge.svg)](http://www.javadoc.io/doc/io.github.andreas-solti.matrix-toolkits-java/mtj)

matrix-toolkits-java 
====================

**MTJ** is a high-performance library for developing linear algebra applications.

See [matrix-toolkits-java](https://github.com/fommil/matrix-toolkits-java) for the original library.

This fork is just a minor extension to cover the case of *general* matrices.
The symmetric case was already covered in [ArpackSym](src/main/java/no/uib/cipr/matrix/sparse/ArpackSym.java).
The new class is [ArpackGen](src/main/java/no/uib/cipr/matrix/sparse/ArpackGen.java).
It uses [ARPACK](http://www.caam.rice.edu/software/ARPACK/)'s [dnaupd](http://www.caam.rice.edu/software/ARPACK/UG/node137.html) and
dneupd routines for the Implicitly Restarted Arnoldi Iteration.


Sparse Solvers
==============

MTJ provides [ARPACK](http://www.caam.rice.edu/software/ARPACK/) for very large symmetric matrices in [ArpackSym](src/main/java/no/uib/cipr/matrix/sparse/ArpackSym.java) (see the example usage in [ArpackSymTest](src/test/java/no/uib/cipr/matrix/sparse/ArpackSymTest.java)). ARPACK solves an arbitrary number of eigenvalues / eigenvectors.

In addition, implementations of the netlib Templates are available in the [`no.uib.cipr.matrix.sparse`](src/test/java/no/uib/cipr/matrix/sparse) package.

Users may wish to look at [Sparse Eigensolvers for Java](http://code.google.com/p/sparse-eigensolvers-java/) for another solver.


Legal
=====

* Copyright (C) 2003-2006 Bjørn-Ove Heimsund
* Copyright (C) 2006-2014 Samuel Halliday


History
=======

This project was originally written by Bjørn-Ove Heimsund, who has taken a step back due to other commitments.
The original project [matrix-toolkits-java](https://github.com/fommil/matrix-toolkits-java) is maintained by Samuel Halliday.


Installation
============

Releases are distributed on Maven central:

```xml
<dependency>
    <groupId>io.github.andreas-solti.matrix-toolkits-java</groupId>
    <artifactId>mtj</artifactId>
    <version>1.0.5</version>
</dependency>
```

Example Code
============
```java
   // check out test class in SparseEigenvalueTest:
   CompColMatrix m = createRandomMatrix(10,15); // create a random
   ArpackGen generalSolver = new ArpackGen(matrix);
   generalSolver.setComputeOnlyEigenvalues(true);
   Map<Double, DenseVectorSub> eigenValueMap = generalSolver.solve(3, ArpackGen.Ritz.LR); // get 3 largest eigenvalues
   double largestARPACKEigenValue = eigenValueMap.keySet().iterator().next();
```


Contributing
============

Contributors are encouraged to fork this repository and issue pull
requests. Contributors implicitly agree to assign an unrestricted licence
to Sam Halliday, but retain the copyright of their code (this means
we both have the freedom to update the licence for those contributions).
