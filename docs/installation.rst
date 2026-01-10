Installation
============

Requirements
------------

asimov-pastro requires Python 3.9 or later.

Core Dependencies
^^^^^^^^^^^^^^^^^

* **asimov** ≥0.5 - Workflow management system
* **numpy** ≥1.20 - Numerical computing
* **scipy** ≥1.7 - Scientific computing
* **h5py** ≥3.0 - HDF5 file support
* **astropy** - Constants and time utilities
* **click** ≥8.0 - Command-line interface

Optional Dependencies
^^^^^^^^^^^^^^^^^^^^^

* **lalsuite** - For accurate antenna pattern calculations (highly recommended)
* **gwpy** ≥3.0 - Gravitational wave data utilities
* **pytest** ≥7.0 - For running tests

Installation Methods
--------------------

From Source (Development)
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   git clone https://git.ligo.org/asimov/pipelines/asimov-pastro.git
   cd asimov-pastro
   pip install -e .

With Development Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   pip install -e ".[dev]"

Verification
------------

Verify the installation:

.. code-block:: bash

   python verify_setup.py

You should see:

.. code-block:: text

   ✓ PASS: Imports
   ✓ PASS: Calculator Creation
   ✓ PASS: Dependencies

   ✓ All tests passed! asimov-pastro is ready to use.

Troubleshooting
---------------

LAL Not Available
^^^^^^^^^^^^^^^^^

If you see warnings about LAL not being available, the pipeline will use simplified antenna pattern calculations. For production use, install LAL:

.. code-block:: bash

   conda install -c conda-forge lalsuite
   # or
   pip install lalsuite

Circular Import Errors
^^^^^^^^^^^^^^^^^^^^^^

If you encounter circular import errors when using the CLI, use the standalone script instead:

.. code-block:: bash

   ./pastro-coherence sample_file.h5

This script loads modules directly without triggering asimov imports.

HDF5 File Reading Issues
^^^^^^^^^^^^^^^^^^^^^^^^^

Ensure you have compatible versions of h5py:

.. code-block:: bash

   pip install --upgrade h5py

For PESummary files, you may need:

.. code-block:: bash

   pip install pesummary
