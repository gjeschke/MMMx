Python
======================


MMMx could be implemented in Matlab, in Python, or as a mixture of both.

Here are some advantages and disadvantages of using Python:

	- Advantages

			- larger potential programmer base
			- larger base of existing bioinformatics modules
			- clean syntax
			- automated documentation is easier
			- fully open source

	- Disadvantages

			- data types are not as good as in Matlab
			- documentation is very heterogeneous across packages 
			- hard to keep consistent, while packages evolve
			- lower programming productivity
			- much overhead bt discussions, which packages should or should not be used
			- porting existing MMM functions is not generally trivial (data typing) 
			- once, everybody programmed Java, now everybody programs Python, will it be Julia tomorrow?


-----------------------


Remarks
-----------------------

* After looking more closely into the language, I am astonished what people can do with it despite its design.

* The syntax is nice and reads well.

* With ``numpy``, Python is well suited to numerical data processing.

* String processing is much better supported in Matlab since 2018

* Compilers have gone out of business, a non-backwards compatible change of the language has been done, there are plans for further changes

* quality and generality of existing modules is very heterogeneous

* many larger projects tend to use Python only as a layer above a C++ core


Current conclusion
------------------

MMMx should be able to integrate methods programmed in Python. These methods can be supplied with *entities*, *rotamer libraries*, *distance distribution restraints*, and *restraint files* by MMMx.

A full port of MMM to Python would delay the first version of MMMx by years.
  