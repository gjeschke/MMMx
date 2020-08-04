.. _property_functions:

Property functions
==========================

Property functions compute properties of an entity or of selected objects within an entity from attribute arrays.

The required attribute arrays can be retrieved by :ref:`object get functions<object_access>`.

Mass
--------------

.. code-block:: matlab

    mass = mass(elements)
	
Parameters
    *   ``elements`` - vector of atomic numbers, (*N*,1) int8 array
Returns
    *   ``mass`` - total mass for the elements vector assuming natural isotope abundance
	 
Use ``get_coor`` with flag ``paradigm`` set to ``true`` for retrieving the element vector.
Alternatively, you can use an address or selection on atom level together with ``get_location`` (slower).