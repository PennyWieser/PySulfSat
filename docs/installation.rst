========================
Installation & Updating
========================

Installation
============

First, obtain Python3 (tested on V3.7). If you haven't used python before, we recomend installing it through anaconda.
 `anaconda3 <https://www.anaconda.com/products/individual>`_.

All of the You Tube examples shown here use Jupyter Lab, which is an easy-friendly editor that will be installed along with Anaconda.

Thermobar can be installed using pip in one line. If you are using a terminal, enter:

.. code-block:: python

   pip install PySCSS

If you are using Jupyter Notebooks or Jupyter Lab, you can also install it by entering the following code into a notebook cell (note the !):

.. code-block:: python

   !pip install PySCSS

You then need to import PySCSS into the script you are running code in. In all the examples, we import Themobar as pt.:

.. code-block:: python

   import PySCSS as ss

This means any time you want to call a function from PySCSS, you do ss.function_name.



Updating
========

To upgrade to the most recent version of PySCSS, type the following into terminal:

.. code-block:: python

   pip install PySCSS --upgrade

Or in your Jupyter environment:

.. code-block:: python

   !pip install Thermobar --upgrade


For maximum reproducability, you should state which version of PySCSS you are using. If you have imported PySCSS as pt, you can find this using:

.. code-block:: python

    ss.__version__