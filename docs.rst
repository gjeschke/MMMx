Documentation
======================

Concept
--------

MMMx documentation will be handled with Sphinx. The main documentation channel is HTML, i.e., the ``.rst`` source files will be optimized for HTML appearance. Still, a PDF version via Latex will be made avialable.  

The documentation is organized as follows:

	- User part
		Everything needed for users without scripting and programming ambitions:

			- introduction to the purpose of MMMx and its structure representation
			- graphical user interface (MMMx Control GUI)
			- command reference

	- Programmer part
		Everything needed for programmatic access, external access, and contributing to MMMx: 

			- introduction to scripting
			- function reference
			- external interface (if developed)
			
An external interface would be implemented as a `REST API'__.

.. __: https://ch.mathworks.com/help/thingspeak/rest-api.html 

Automation
--------------------------------

The following aspects will be automated in order to keep MMMx documentation and package consistent:

* basic description of GUI windows

* command reference

* function reference

For these cases, the documentation source are comment sections in the code source files. Compilation of the documentation has two steps:

1) extraction of the information from the code source files into ``.rst`` source files

2) building the HTML (and Latex) files with the ``make`` command of Sphinx

For Python code, the first step is done by the ``autodoc`` feature of Sphinx. For Matlab code, this will be done by a Matlab function.


