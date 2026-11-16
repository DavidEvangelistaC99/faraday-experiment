# Signal Chain

Signal Chain is a radar data processing library wich includes modules to read,
and write different files formats, besides modules to process and visualize the
data.

## Dependencies

- GCC (gcc and gfortran)
- Python.h (python-dev or python-devel)
- Python-TK (python-tk)
- HDF5 libraries (libhdf5-dev)

## Installation

### From source

First, ensure that you have the above-listed dependencies installed, then clone
the repository and following next steps:

```
conda create -n schain-proc python=3.6
conda activate schain-proc
pip install pip==23.3.1
pip install "setuptools<60"
pip install "numpy<1.24"
pip install -e .
```

Alternative:

```
python setup.py develop
```

Next:

```
pip install "click<8.2"
pip install "numpy<1.24"
pip install zmq
pip install matplotlib
pip install scipy
pip install digital_rf==2.6.7
pip install cartopy
pip install mysql-connector-python
pip install pandas
```

## Example

Here you can find an script to read Spectra data (.pdata), remove dc and plot
self-spectra & RTI:

```python
#!/usr/bin/python

from schainpy.controller import Project

prj = Project()

read_unit = prj.addReadUnit(
    datatype='Spectra',
    path='/path/to/pdata/',
    startDate='2014/01/31',
    endDate='2014/03/31',
    startTime='00:00:00',
    endTime='23:59:59',
    online=0,
    walk=0
    )

proc_unit = prj.addProcUnit(
    datatype='Spectra',
    inputId=read_unit.getId()
    )

op = proc_unit.addOperation(name='selectChannels')
op.addParameter(name='channelList', value='0,1')

op = proc_unit.addOperation(name='selectHeights')
op.addParameter(name='minHei', value='80')
op.addParameter(name='maxHei', value='200')

op = proc_unit.addOperation(name='removeDC')

op = proc_unit.addOperation(name='SpectraPlot')
op.addParameter(name='wintitle', value='Spectra', format='str')

op = proc_unit.addOperation(name='RTIPlot')
op.addParameter(name='wintitle', value='RTI', format='str')

prj.start()

```
